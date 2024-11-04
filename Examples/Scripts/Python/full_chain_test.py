#!/usr/bin/env python3

import sys, os, argparse, pathlib
import acts, acts.examples


def parse_args():
    from acts.examples.reconstruction import SeedingAlgorithm

    parser = argparse.ArgumentParser(
        description="Script to test full chain ACTS simulation and reconstruction",
    )
    parser.add_argument(
        "-G",
        "--generic-detector",
        action="store_true",
        help="Use generic detector geometry and config",
    )
    parser.add_argument(
        "--odd",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Use Open Data Detector geometry and config (default unless overridden by -G or -A). Requires ACTS_BUILD_ODD.",
    )
    parser.add_argument(
        "-A",
        "--itk",
        action="store_true",
        help="Use ATLAS ITk geometry and config. Requires acts-itk/ in current directory.",
    )
    parser.add_argument(
        "-n",
        "--events",
        type=int,
        default=100,
        help="The number of events to process (default=%(default)d).",
    )
    parser.add_argument(
        "-s",
        "--skip",
        type=int,
        default=0,
        help="Number of events to skip (default=%(default)d)",
    )
    parser.add_argument(
        "-N",
        "--gen-nparticles",
        type=int,
        default=4,
        help="Number of generated particles per vertex from the particle gun (default=%(default)d).",
    )
    parser.add_argument(
        "-M",
        "--gen-nvertices",
        type=int,
        default=200,
        help="Number of vertices per event (multiplicity) from the particle gun; or number of pileup events (default=%(default)d)",
    )
    parser.add_argument(
        "-j",
        "--threads",
        type=int,
        default=-1,
        help="Number of parallel threads, negative for automatic (default).",
    )
    parser.add_argument(
        "-t",
        "--ttbar-pu200",
        action="store_true",
        help="Generate ttbar + mu=200 pile-up using Pythia8",
    )
    parser.add_argument(
        "-r",
        "--random-seed",
        type=int,
        default=42,
        help="Random number seed (default=%(default)d)",
    )
    parser.add_argument(
        "-l",
        "--loglevel",
        type=int,
        default=2,
        help="The output log level. Please set the wished number (0 = VERBOSE, 1 = DEBUG, 2 = INFO (default), 3 = WARNING, 4 = ERROR, 5 = FATAL).",
    )
    parser.add_argument(
        "-p",
        "--gen-mom-gev",
        default="1:10",
        help="pT - transverse momentum generation range in GeV (min:max pT for particle gun; min pT for -t) (default 1:10, except max Pt not used with -t)",
    )
    parser.add_argument(
        "--digi-config",
        type=pathlib.Path,
        help="Digitization configuration file",
    )
    parser.add_argument(
        "--material-config",
        type=pathlib.Path,
        help="Material map configuration file",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        "--output",
        default=None,
        type=pathlib.Path,
        help="Directory to write outputs to",
    )
    parser.add_argument(
        "-O",
        "--output-detail",
        action="count",
        default=0,
        help="fewer output files. Use -OO for more output files. Use -OOO to disable all output.",
    )
    parser.add_argument(
        "-d",
        "--dump-args-calls",
        action="store_true",
        help="Show pybind function call details",
    )
    parser.add_argument(
        "-S",
        "--seeding-algorithm",
        action=EnumAction,
        enum=SeedingAlgorithm,
        default=SeedingAlgorithm.Default,
        help="Select the seeding algorithm to use",
    )
    parser.add_argument(
        "-e",
        "--gen-cos-theta",
        action="store_true",
        help="Sample eta as cos(theta) and not uniform",
    )
    parser.add_argument(
        "-b",
        "--bf-constant",
        action="store_true",
        help="Use constant 2T B-field also for ITk; and don't include material map",
    )
    parser.add_argument(
        "--edm4hep",
        type=pathlib.Path,
        help="Use edm4hep inputs",
    )
    parser.add_argument(
        "-g",
        "--geant4",
        action="store_true",
        help="Use Geant4 instead of Fatras for detector simulation",
    )
    parser.add_argument(
        "-F",
        "--disable-fpemon",
        action="store_true",
        help="sets ACTS_SEQUENCER_DISABLE_FPEMON=1",
    )
    parser.add_argument(
        "--ckf",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Switch CKF on/off",
    )
    parser.add_argument(
        "--reco",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Switch reco on/off",
    )
    parser.add_argument(
        "--vertexing",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Switch vertexing on/off",
    )
    parser.add_argument(
        "--MLSeedFilter",
        action="store_true",
        help="Use the Ml seed filter to select seed after the seeding step",
    )
    parser.add_argument(
        "-C",
        "--simple-ckf",
        action="store_true",
        help="Turn off CKF features: seed deduplication, two-way CKF, stick on the seed measurements during track finding, max pixel/strip holes",
    )
    parser.add_argument(
        "--ambi-solver",
        type=str,
        choices=["greedy", "scoring", "ML", "none"],
        default="greedy",
        help="Set which ambiguity solver to use (default=%(default)s)",
    )
    parser.add_argument(
        "--ambi-config",
        type=pathlib.Path,
        default=pathlib.Path.cwd() / "ambi_config.json",
        help="Set the configuration file for the Score Based ambiguity resolution",
    )
    parser.add_argument(
        "-c",
        "--output-csv",
        action="count",
        default=0,
        help="Use CSV output instead of ROOT. Specify -cc to output both.",
    )
    return parser.parse_args()


def full_chain(args):
    import acts

    # keep these in memory after we return the sequence
    global detector, trackingGeometry, decorators, field, rnd

    if args.dump_args_calls:
        acts.examples.dump_args_calls(locals())

    logger = acts.logging.getLogger("full_chain_test")

    nDetArgs = [args.generic_detector, args.odd, args.itk].count(True)
    if nDetArgs == 0:
        args.generic_detector = True
    elif nDetArgs == 2:
        args.odd = False
    nDetArgs = [args.generic_detector, args.odd, args.itk].count(True)
    if nDetArgs != 1:
        logger.fatal("require exactly one of: --generic-detector --odd --itk")
        sys.exit(2)
    if args.generic_detector:
        detname = "gen"
    elif args.itk:
        detname = "itk"
    elif args.odd:
        detname = "odd"

    u = acts.UnitConstants
    pt = [float(p) * u.GeV if p != "" else None for p in args.gen_mom_gev.split(":")]
    if len(pt) == 1:
        pt.append(pt[0])  # 100 -> 100:100
    if len(pt) != 2:
        logger.fatal(f"bad option value: --gen-mom-gev {args.gen_mom_gev}")
        sys.exit(2)

    if args.output_detail == 3:
        outputDirLess = None
    elif args.output_dir is None:
        outputDirLess = pathlib.Path.cwd() / f"{detname}_output"
    else:
        outputDirLess = args.output_dir

    outputDir = None if args.output_detail == 1 else outputDirLess
    outputDirMore = None if args.output_detail in (0, 1) else outputDirLess

    outputDirRoot = outputDir if args.output_csv != 1 else None
    outputDirLessRoot = outputDirLess if args.output_csv != 1 else None
    outputDirMoreRoot = outputDirMore if args.output_csv != 1 else None
    outputDirCsv = outputDir if args.output_csv != 0 else None
    outputDirLessCsv = outputDirLess if args.output_csv != 0 else None
    outputDirMoreCsv = outputDirMore if args.output_csv != 0 else None

    if args.disable_fpemon:
        os.environ["ACTS_SEQUENCER_DISABLE_FPEMON"] = "1"

    # fmt: off
    if args.generic_detector:
        etaRange = (-2.0, 2.0)
        rhoMax = 24.0 * u.mm
        geo_dir = pathlib.Path(acts.__file__).resolve().parent.parent.parent.parent.parent
        if args.loglevel <= 2: logger.info(f"Load Generic Detector from {geo_dir}")
        if args.digi_config is None:
            args.digi_config = geo_dir / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json"
        seedingConfigFile = geo_dir / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
        args.bf_constant = True
        detector, trackingGeometry, decorators = acts.examples.GenericDetector.create()
    elif args.odd:
        import acts.examples.odd
        etaRange = (-3.0, 3.0)
        rhoMax = 24.0 * u.mm
        beamTime = 1.0 * u.ns
        geo_dir = acts.examples.odd.getOpenDataDetectorDirectory()
        if args.loglevel <= 2: logger.info(f"Load Open Data Detector from {geo_dir.resolve()}")
        if args.digi_config is None:
            args.digi_config = geo_dir / "config/odd-digi-smearing-config.json"
        seedingConfigFile = geo_dir / "config/odd-seeding-config.json"
        if args.material_config is None:
            args.material_config = geo_dir / "data/odd-material-maps.root"
        args.bf_constant = True
        detector, trackingGeometry, decorators = acts.examples.odd.getOpenDataDetector(
            odd_dir=geo_dir,
            mdecorator=acts.IMaterialDecorator.fromFile(args.material_config)
            if not args.bf_constant
            else None,
        )
    elif args.itk:
        import acts.examples.itk as itk
        etaRange = (-4.0, 4.0)
        rhoMax = 28.0 * u.mm
        beamTime = 5.0 * u.ns
        geo_dir = pathlib.Path("acts-itk")
        if args.loglevel <= 2: logger.info(f"Load ATLAS ITk from {geo_dir.resolve()}")
        if args.digi_config is None:
            args.digi_config = geo_dir / "itk-hgtd/itk-smearing-config.json"
        seedingConfigFile = geo_dir / "itk-hgtd/geoSelection-ITk.json"
        # args.material_config defaulted in itk.buildITkGeometry: geo_dir / "itk-hgtd/material-maps-ITk-HGTD.json"
        bFieldFile = geo_dir / "bfield/ATLAS-BField-xyz.root"
        detector, trackingGeometry, decorators = itk.buildITkGeometry(
            geo_dir,
            customMaterialFile=args.material_config,
            material=not args.bf_constant,
            logLevel=acts.logging.Level(args.loglevel),
        )
    # fmt: on

    if args.bf_constant:
        field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
    else:
        logger.info("Create magnetic field map from %s" % str(bFieldFile))
        field = acts.examples.MagneticFieldMapXyz(str(bFieldFile))
    rnd = acts.examples.RandomNumbers(seed=42)

    from acts.examples.simulation import (
        MomentumConfig,
        EtaConfig,
        ParticleConfig,
        ParticleSelectorConfig,
        addDigitization,
        addParticleSelection,
    )

    s = acts.examples.Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=args.threads if not (args.geant4 and args.threads == -1) else 1,
        logLevel=acts.logging.Level(args.loglevel),
        outputDir="" if outputDir is None else str(outputDir),
    )

    # is this needed?
    for d in decorators:
        s.addContextDecorator(d)

    preSelectParticles = (
        ParticleSelectorConfig(
            rho=(0.0 * u.mm, rhoMax),
            absZ=(0.0 * u.mm, 1.0 * u.m),
            eta=etaRange,
            pt=(150 * u.MeV, None),
            removeNeutral=True,
        )
        if args.edm4hep or args.geant4 or args.ttbar_pu200
        else ParticleSelectorConfig()
    )

    postSelectParticles = (
        ParticleSelectorConfig(
            pt=(pt[0], None),
            eta=etaRange,
            measurements=(9, None),
            removeNeutral=True,
        )
        if args.ttbar_pu200
        else ParticleSelectorConfig()
    )

    if args.edm4hep:
        import acts.examples.edm4hep

        edm4hepReader = acts.examples.edm4hep.EDM4hepReader(
            inputPath=str(args.edm4hep),
            inputSimHits=[
                "PixelBarrelReadout",
                "PixelEndcapReadout",
                "ShortStripBarrelReadout",
                "ShortStripEndcapReadout",
                "LongStripBarrelReadout",
                "LongStripEndcapReadout",
            ],
            outputParticlesGenerator="particles_input",
            outputParticlesInitial="particles_initial",
            outputParticlesFinal="particles_final",
            outputSimHits="simhits",
            graphvizOutput="graphviz",
            dd4hepDetector=detector,
            trackingGeometry=trackingGeometry,
            sortSimHitsInTime=True,
            level=acts.logging.INFO,
        )
        s.addReader(edm4hepReader)
        s.addWhiteboardAlias("particles", edm4hepReader.config.outputParticlesGenerator)

        addParticleSelection(
            s,
            config=preSelectParticles,
            inputParticles="particles",
            outputParticles="particles_selected",
        )

    else:

        if not args.ttbar_pu200:
            from acts.examples.simulation import addParticleGun

            addParticleGun(
                s,
                MomentumConfig(*pt, transverse=True),
                EtaConfig(*etaRange, uniform=not args.gen_cos_theta),
                ParticleConfig(
                    args.gen_nparticles, acts.PdgParticle.eMuon, randomizeCharge=True
                ),
                multiplicity=args.gen_nvertices,
                rnd=rnd,
                outputDirRoot=outputDirMoreRoot,
                outputDirCsv=outputDirMoreCsv,
            )
        else:
            from acts.examples.simulation import addPythia8

            addPythia8(
                s,
                hardProcess=["Top:qqbar2ttbar=on"],
                npileup=args.gen_nvertices,
                vtxGen=acts.examples.GaussianVertexGenerator(
                    stddev=acts.Vector4(
                        0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                    ),
                    mean=acts.Vector4(0, 0, 0, 0),
                ),
                rnd=rnd,
                outputDirRoot=outputDirLessRoot,
                outputDirCsv=outputDirLessCsv,
            )

        if not args.geant4:
            from acts.examples.simulation import addFatras

            addFatras(
                s,
                trackingGeometry,
                field,
                rnd=rnd,
                preSelectParticles=preSelectParticles,
                postSelectParticles=postSelectParticles,
                outputDirRoot=outputDirLessRoot,
                outputDirCsv=outputDirLessCsv,
            )
        else:
            if s.config.numThreads != 1:
                logger.fatal(
                    f"Geant 4 simulation does not support multi-threading (threads={s.config.numThreads})"
                )
                sys.exit(2)

            from acts.examples.simulation import addGeant4

            # Pythia can sometime simulate particles outside the world volume, a cut on the Z of the track help mitigate this effect
            # Older version of G4 might not work, this as has been tested on version `geant4-11-00-patch-03`
            # For more detail see issue #1578
            addGeant4(
                s,
                detector,
                trackingGeometry,
                field,
                rnd=rnd,
                preSelectParticles=preSelectParticles,
                postSelectParticles=postSelectParticles,
                killVolume=trackingGeometry.highestTrackingVolume,
                killAfterTime=25 * u.ns,
                outputDirRoot=outputDirLessRoot,
                outputDirCsv=outputDirLessCsv,
            )

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=args.digi_config,
        rnd=rnd,
        outputDirRoot=outputDirLessRoot,
        outputDirCsv=outputDirLessCsv,
    )

    if not args.reco:
        return s

    from acts.examples.reconstruction import (
        addSeeding,
        ParticleSmearingSigmas,
        addCKFTracks,
        CkfConfig,
        SeedingAlgorithm,
        TrackSelectorConfig,
        addAmbiguityResolution,
        AmbiguityResolutionConfig,
        addVertexFitting,
        VertexFinder,
    )

    if args.itk and args.seeding_algorithm == SeedingAlgorithm.Default:
        seedingAlgConfig = itk.itkSeedingAlgConfig(
            itk.InputSpacePointsType.PixelSpacePoints
        )
    else:
        seedingAlgConfig = []

    addSeeding(
        s,
        trackingGeometry,
        field,
        *seedingAlgConfig,
        ParticleSmearingSigmas(
            ptRel=0.01
        ),  # only needed for SeedingAlgorithm.TruthSmeared
        seedingAlgorithm=args.seeding_algorithm,
        rnd=rnd,  # only needed for SeedingAlgorithm.TruthSmeared
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0.1 * u.e / u.GeV,
            1 * u.ns,
        ],
        initialSigmaPtRel=0.1,
        initialVarInflation=[1.0] * 6,
        geoSelectionConfigFile=seedingConfigFile,
        outputDirRoot=outputDirRoot,
        outputDirCsv=outputDirCsv,
    )

    if args.MLSeedFilter:
        from acts.examples.reconstruction import (
            addSeedFilterML,
            SeedFilterMLDBScanConfig,
        )

        addSeedFilterML(
            s,
            SeedFilterMLDBScanConfig(
                epsilonDBScan=0.03, minPointsDBScan=2, minSeedScore=0.1
            ),
            onnxModelFile=str(
                geo_dir
                / "Examples/Scripts/Python/MLAmbiguityResolution/seedDuplicateClassifier.onnx"
            ),
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
        )

    if not args.ckf:
        return s

    if not args.itk:

        trackSelectorConfig = TrackSelectorConfig(
            pt=(500 * u.MeV, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
            nMeasurementsMin=7,
            maxHoles=2,
            maxOutliers=2,
        )
        # fmt: off
        ckfConfig = CkfConfig(
            chi2CutOffMeasurement=15.0,
            chi2CutOffOutlier=25.0,
            numMeasurementsCutOff=10,
            **(dict(
                seedDeduplication=True,
                stayOnSeed=True,
            ) if not args.simple_ckf and args.seeding_algorithm != SeedingAlgorithm.TruthSmeared else {}),
            **(dict(
                pixelVolumes=[16, 17, 18],
                stripVolumes=[23, 24, 25],
                maxPixelHoles=1,
                maxStripHoles=2,
                constrainToVolumes=[
                    2,  # beam pipe
                    32,
                    4,  # beam pip gap
                    16,
                    17,
                    18,  # pixel
                    20,  # PST
                    23,
                    24,
                    25,  # short strip
                    26,
                    8,  # long strip gap
                    28,
                    29,
                    30,  # long strip
                ],
            ) if args.odd and not args.simple_ckf else {})
        )
        # fmt: on

    elif args.itk:
        # fmt: off
        trackSelectorConfig = (
            TrackSelectorConfig(absEta=(None, 2.0), pt=(0.9 * u.GeV, None), nMeasurementsMin=9, maxHoles=2, maxOutliers=2, maxSharedHits=2),
            TrackSelectorConfig(absEta=(None, 2.6), pt=(0.4 * u.GeV, None), nMeasurementsMin=8, maxHoles=2, maxOutliers=2, maxSharedHits=2),
            TrackSelectorConfig(absEta=(None, 4.0), pt=(0.4 * u.GeV, None), nMeasurementsMin=7, maxHoles=2, maxOutliers=2, maxSharedHits=2),
        )
        ckfConfig = CkfConfig(
            **(dict(
                seedDeduplication=True,
                stayOnSeed=True,
            ) if args.seeding_algorithm != SeedingAlgorithm.TruthSmeared else {}),
            # ITk volumes from Noemi's plot
            pixelVolumes=[8, 9, 10, 13, 14, 15, 16, 18, 19, 20],
            stripVolumes=[22, 23, 24],
            maxPixelHoles=1,
            maxStripHoles=2,
        ) if not args.simple_ckf else CkfConfig()
        # fmt: on

    if args.output_detail == 1:
        writeDetail = dict(writeTrackSummary=False)
    elif args.output_detail == 2:
        writeDetail = dict(writeTrackStates=True)
    else:
        writeDetail = {}

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        trackSelectorConfig=trackSelectorConfig,
        ckfConfig=ckfConfig,
        **(dict(twoWay=False) if not args.simple_ckf else {}),
        **writeDetail,
        outputDirRoot=outputDirRoot,
        outputDirCsv=outputDirCsv,
    )

    if args.ambi_solver == "ML":

        from acts.examples.reconstruction import (
            addAmbiguityResolutionML,
            AmbiguityResolutionMLConfig,
        )

        addAmbiguityResolutionML(
            s,
            AmbiguityResolutionMLConfig(
                maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
            ),
            onnxModelFile=str(
                geo_dir
                / "Examples/Scripts/Python/MLAmbiguityResolution/duplicateClassifier.onnx"
            ),
            outputDirRoot=outputDirLessRoot,
            outputDirCsv=outputDirLessCsv,
        )

    elif args.ambi_solver == "scoring":

        from acts.examples.reconstruction import (
            addScoreBasedAmbiguityResolution,
            ScoreBasedAmbiguityResolutionConfig,
        )
        import math

        addScoreBasedAmbiguityResolution(
            s,
            ScoreBasedAmbiguityResolutionConfig(
                minScore=0,
                minScoreSharedTracks=1,
                maxShared=2,
                maxSharedTracksPerMeasurement=2,
                pTMax=1400,
                pTMin=0.5,
                phiMax=math.pi,
                phiMin=-math.pi,
                etaMax=4,
                etaMin=-4,
                useAmbiguityFunction=False,
            ),
            ambiVolumeFile=args.ambi_config,
            **(dict(writeCovMat=True) if args.output_detail != 1 else {}),
            outputDirRoot=outputDirLessRoot,
            outputDirCsv=outputDirLessCsv,
        )

    elif args.ambi_solver == "greedy":

        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=3,
                maximumIterations=10000 if args.itk else 1000000,
                nMeasurementsMin=6 if args.itk else 7,
            ),
            **writeDetail,
            outputDirRoot=outputDirLessRoot,
            outputDirCsv=outputDirLessCsv,
        )

    if args.vertexing:
        addVertexFitting(
            s,
            field,
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=outputDirLess,
        )

    return s


# Graciously taken from https://stackoverflow.com/a/60750535/4280680 (via seeding.py)
class EnumAction(argparse.Action):
    """
    Argparse action for handling Enums
    """

    def __init__(self, **kwargs):
        import enum

        # Pop off the type value
        enum_type = kwargs.pop("enum", None)

        # Ensure an Enum subclass is provided
        if enum_type is None:
            raise ValueError("type must be assigned an Enum when using EnumAction")
        if not issubclass(enum_type, enum.Enum):
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


# use modules to prevent Fluke unused import warning
dir(acts)
dir(acts.examples)

# main program: parse arguments, setup sequence, and run the full chain
full_chain(parse_args()).run()
