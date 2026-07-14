import acts
import argparse
from acts.examples.simulation import (
    addParticleGun,
    addGeant4,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    MomentumConfig,
)

from acts.examples.odd import getOpenDataDetector

u = acts.UnitConstants


def main():

    import acts

    global __geant4Handle
    global __pmr

    p = argparse.ArgumentParser()

    p.add_argument("-n", "--events", type=int, default=1000, help="Number of Events")

    p.add_argument(
        "-j",
        "--threads",
        type=int,
        default=1,
        help="Number of Threads for parallel execution",
    )

    p.add_argument("-i", "--input", type=str, default="", help="Input SQL file")

    p.add_argument("-o", "--output", type=str, default="", help="Output prefix")

    p.add_argument(
        "-m", "--map", type=str, default="", help="Input file for the material map"
    )

    p.add_argument("-s", "--seed", type=int, default=221177, help="Random number seed")

    # Add Particle generation related arguments
    p.add_argument(
        "--tracks", type=int, default=10, help="Number of tracks to generate per event"
    )
    p.add_argument(
        "--eta-range",
        type=float,
        nargs=2,
        default=[-2.5, 2.5],
        help="Pseudorapidity range for particle generation",
    )
    p.add_argument(
        "--phi-range",
        type=float,
        nargs=2,
        default=[0, 360],
        help="Azimuthal angle range (in degrees) for particle generation",
    )
    p.add_argument(
        "--pt-range",
        type=float,
        nargs=2,
        default=[0.5, 10],
        help="Transverse momentum range (in GeV) for particle generation",
    )

    # The geometry modes are
    # gen1: Gen1 detector with Gen1 navigator and propagator
    # gen3: Gen3 detector with Gen3 navigator and propagator
    # detray: Gen3->detray detector with detray navigator and propagator
    # geant4: Geant4 navigator and propagator with gen3 surface matching
    p.add_argument(
        "--geo-mode",
        type=str,
        default="gen3",
        choices=["gen1", "gen3", "detray", "geant4"],
        help="Convert to detray detector and run detray navigation and propagation",
    )

    # Output options

    p.add_argument(
        "--output-summary",
        action=argparse.BooleanOptionalAction,
        help="Write out the summary objects",
    )

    p.add_argument(
        "--output-steps",
        action=argparse.BooleanOptionalAction,
        help="Write out the step objects",
    )

    p.add_argument(
        "--output-material",
        action=argparse.BooleanOptionalAction,
        help="Write out the recorded",
    )

    p.add_argument(
        "--output-sim-hits",
        action=argparse.BooleanOptionalAction,
        help="Write out sim hits, only makes sense for Geant4",
    )

    p.add_argument(
        "--detray-search-window",
        type=int,
        nargs=2,
        default=[0, 0],
        metavar=("AXIS0", "AXIS1"),
        help="Detray grid acceleration search window per axis. {0, 0} (default) "
        "only looks at the current bin; increase to probe whether the converted "
        "surface grids are missing neighbor-cell entries",
    )

    args = p.parse_args()

    prfx = args.output + "_" if args.output != "" else ""

    gContext = acts.GeometryContext.dangerouslyDefaultConstruct()
    logLevel = acts.logging.INFO

    # Common (to all modes): Evoke the sequence
    rnd = acts.examples.RandomNumbers(seed=args.seed)

    # Common: Build the sequencer
    s = acts.examples.Sequencer(events=args.events, numThreads=args.threads)

    # Common: Add the particle gun from ACTS
    addParticleGun(
        s,
        ParticleConfig(
            num=args.tracks, pdg=acts.PdgParticle.eMuon, randomizeCharge=True
        ),
        EtaConfig(args.eta_range[0], args.eta_range[1]),
        PhiConfig(args.phi_range[0], args.phi_range[1]),
        MomentumConfig(
            args.pt_range[0] * u.GeV, args.pt_range[1] * u.GeV, transverse=True
        ),
        rnd=rnd,
    )

    # Timing measurement is run if neither output in on
    sterileRun = False
    if not args.output_summary and not args.output_steps and not args.output_material:
        print(">> Timing measurement is enabled, no output is written")
        sterileRun = True

    # Material decoration for reconstruction geometry
    materialDecorator = None
    if args.map != "":
        print(">>> Loading a material decorator from file:", args.map)
        materialDecorator = acts.IMaterialDecorator.fromFile(args.map)

    trackingGeometry = None
    detectorStore = {}

    # Build from ODD unless reading detray geometry from file
    build_from_odd = args.input == "" or args.geo_mode != "detray"

    if build_from_odd:
        buildGen3 = args.geo_mode == "gen3" or args.geo_mode == "detray"
        with getOpenDataDetector(gen3=buildGen3) as detector:
            trackingGeometry = detector.trackingGeometry()
            detectorStore["Detector"] = detector
            detectorStore["Volume"] = trackingGeometry.highestTrackingVolume
            detectorStore["SurfaceByIdentifier"] = trackingGeometry.geoIdSurfaceMap()

    print(">>> Test mode is :", args.geo_mode)
    # check if the mode does not contain geant4
    if args.geo_mode != "geant4":
        # The propagator
        propagatorImpl = None
        stepper = acts.StraightLineStepper()

        # Build the detector for Gen1/Gen3
        if args.geo_mode != "detray":
            # Set up the navigator - Gen1/Gen3
            navigator = acts.Navigator(trackingGeometry=trackingGeometry)
            propagator = acts.Propagator(stepper, navigator)
            propagatorImpl = acts.examples.ConcretePropagator(propagator)
        else:
            import glob
            import acts.vecmem, acts.detray
            import acts.examples.detray

            __pmr = acts.vecmem.HostMemoryResource()

            if args.input != "":
                files = glob.glob(args.input.rstrip("/") + "/*.json")
                print(">>> Reading detray geometry from", args.input, "->", files)
                detrayGeometry, _ = acts.detray.readODD(__pmr, files)
            else:
                detrayGeometry, _ = acts.detray.convertODD(
                    __pmr,
                    gContext,
                    trackingGeometry,
                    beampipeVolumeName="BeamPipe",
                    logLevel=logLevel,
                    convertMaterial=True,
                    convertSurfaceGrids=True,
                )
            propagatorImpl = acts.examples.detray.StraightLinePropagatorODD(
                detrayGeometry,
                __pmr,
                sterileRun,
                logLevel,
                searchWindow=args.detray_search_window,
            )

        # Run particle smearing
        trkParamExtractor = acts.examples.ParticleTrackParamExtractor(
            level=acts.logging.INFO,
            inputParticles="particles_generated",
            outputTrackParameters="start_parameters",
        )
        s.addAlgorithm(trkParamExtractor)

        propagationAlgorithm = acts.examples.PropagationAlgorithm(
            propagatorImpl=propagatorImpl,
            level=acts.logging.INFO,
            sterileLogger=sterileRun,
            inputTrackParameters="start_parameters",
            outputSummaryCollection="propagation_summary",
            outputMaterialCollection="material_tracks",
        )
        s.addAlgorithm(propagationAlgorithm)
    else:
        print(">>> Running Geant4 simulation, buckle up...")
        detector = detectorStore["Detector"]

        import acts.examples.geant4

        from acts.examples.geant4 import (
            Geant4Simulation,
            Geant4ConstructionOptions,
            SensitiveSurfaceMapper,
        )

        customLogLevel = acts.examples.defaultLogging(s, logLevel)

        smmConfig = SensitiveSurfaceMapper.Config()
        smmConfig.volumeMappings = []
        smmConfig.materialMappings = ["Silicon"]
        sensitiveMapper = SensitiveSurfaceMapper.create(
            smmConfig, customLogLevel(), trackingGeometry
        )

        detectorConstructionOptions = acts.examples.geant4.Geant4ConstructionOptions()

        # Specify the physics list, kill volume and other options for the Geant4 simulation
        physicsList = "FTFP_BERT"  # "MaterialPhysicsList"
        killVolume = detectorStore["Volume"]
        killAfterTime = float("inf")
        bfield = None

        inputParticles = "particles_generated"
        outputParticles = "particles_final"

        addGeant4(
            s,
            detector,
            trackingGeometry,
            None,
            outputDirRoot=None,
            outputDirCsv=None,
            outputDirObj=None,
            rnd=rnd,
            killVolume=trackingGeometry.highestTrackingVolume,
            killAfterTime=25 * u.ns,
            killSecondaries=True,
            recordHitsOfSecondaries=False,
            recordPropagationSummaries=True,
        )

        s.addWhiteboardAlias("propagation_summary", "propagation_summaries")

    # Common: Write the summary
    if args.output_summary:
        s.addWriter(
            acts.examples.root.RootPropagationSummaryWriter(
                level=acts.logging.INFO,
                inputSummaryCollection="propagation_summary",
                filePath=prfx + args.geo_mode + "_propagation_summary.root",
            )
        )

    # Common: Write the steps
    if args.output_steps:
        s.addWriter(
            acts.examples.root.RootPropagationStepsWriter(
                level=acts.logging.INFO,
                collection="propagation_summary",
                filePath=prfx + args.geo_mode + "_propagation_steps.root",
            )
        )

    # Common: Write the material
    if args.output_material:
        s.addWriter(
            acts.examples.root.RootMaterialTrackWriter(
                level=acts.logging.INFO,
                inputMaterialTracks="material_tracks",
                filePath=args.geo_mode + "_material_tracks.root",
                storeSurface=False,
                storeVolume=False,
            )
        )

    # Run the sequence
    s.run()


if "__main__" == __name__:
    main()
