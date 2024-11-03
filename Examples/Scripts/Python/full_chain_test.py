#!/usr/bin/env python3

import sys, os

# needed if this script is a symlink to another directory
sys.path.insert(0, os.path.dirname(__file__))

import pathlib, argparse, acts, acts.examples, seeding, acts.examples.reconstruction

parser = argparse.ArgumentParser(
    description="Script to test full chain ACTS simulation and reconstruction",
)
parser.add_argument(
    "-n",
    "--events",
    type=int,
    default=100,
    help="The number of events to process (default=100).",
)
parser.add_argument(
    "-s",
    "--skip",
    type=int,
    default=0,
    help="Number of events to skip (default=0)",
)
parser.add_argument(
    "-N",
    "--gen-nparticles",
    type=int,
    help="Number of generated particles per event (default=2) or number of pileup events (default=200).",
)
parser.add_argument(
    "--gen-multiplicity",
    type=int,
    default=1,  # full_chain_odd.py has 200
    help="Multiplicity (no. of vertices) of the particle gun",
)
parser.add_argument(
    "-j",
    "--jobs",
    type=int,
    default=-1,
    help="Number of parallel jobs, negative for automatic (default).",
)
parser.add_argument(
    "-t",
    "--ttbar-pu200",
    action="store_true",
    help="Generate ttbar + mu=200 pile-up using Pythia8",
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
    help="pT - transverse momentum generation range in GeV (default 1:10 or 1: with -t)",
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
    "--output-dir",
    "--output",
    "-o",
    default=None,
    type=pathlib.Path,
    help="Directory to write outputs to",
)
parser.add_argument(
    "-a",
    "--algorithm",
    action=seeding.EnumAction,
    enum=acts.examples.reconstruction.SeedingAlgorithm,
    default=acts.examples.reconstruction.SeedingAlgorithm.Default,
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
    "-G",
    "--generic-detector",
    action="store_true",
    help="Use generic detector geometry and config",
)
parser.add_argument(
    "--odd",
    default=True,
    action=argparse.BooleanOptionalAction,
    help="Use Open Data Detector geometry and config (default unless overridden by -G or -A)",
)
parser.add_argument(
    "-A",
    "--itk",
    action="store_true",
    help="Use ATLAS ITk geometry and config",
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
    "-d",
    "--dump-args-calls",
    action="store_true",
    help="Show pybind function call details",
)
parser.add_argument(
    "-O",
    "--reduce-output",
    action="count",
    default=0,
    help="don't write intermediate results. Use -OO to disable all output.",
)
parser.add_argument(
    "-F",
    "--disable-fpemon",
    action="store_true",
    help="sets ACTS_SEQUENCER_DISABLE_FPEMON=1",
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
    help="Set which ambiguity solver to use (default=greedy)",
)
parser.add_argument(
    "--ambi-config",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "ambi_config.json",
    help="Set the configuration file for the Score Based ambiguity resolution",
)
parser.add_argument(
    "--output-csv",
    action="store_true",
    help="Use CSV output instead of ROOT",
)
args = parser.parse_args()

if args.disable_fpemon:
    os.environ["ACTS_SEQUENCER_DISABLE_FPEMON"] = "1"

u = acts.UnitConstants
if args.itk or not (args.odd or args.generic_detector):
    detname = "itk"
    args.odd = False
    args.generic_detector = False
if args.generic_detector:
    detname = "gen"
    args.odd = False
if args.odd:
    detname = "odd"
if args.gen_mom_gev is None:
    # didn't need to do this, because max pT isn't used for ttbar
    args.gen_mom_gev = "1:" if args.ttbar_pu200 else "1:10"
pt = [float(p) * u.GeV if p != "" else None for p in args.gen_mom_gev.split(":")]
if args.gen_nparticles is None:
    args.gen_nparticles = 200 if args.ttbar_pu200 else 2

if args.dump_args_calls:
    acts.examples.dump_args_calls(locals())

logger = acts.logging.getLogger("full_chain_test")

if args.reduce_output >= 2:
    outputDirFinal = None
elif args.output_dir is None:
    outputDirFinal = pathlib.Path.cwd() / f"{detname}_output"
else:
    outputDirFinal = args.output_dir
outputDir = None if args.reduce_output >= 1 else outputDirFinal

outputDirRoot = outputDir if not args.output_csv else None
outputDirFinalRoot = outputDirFinal if not args.output_csv else None
outputDirCsv = outputDir if args.output_csv else None
outputDirFinalCsv = outputDirFinal if args.output_csv else None

# fmt: off
if args.generic_detector:
    etaRange = (-2.0, 2.0)
    rhoMax = 24.0 * u.mm
    geo_dir = pathlib.Path(acts.__file__).resolve().parent.parent.parent.parent.parent
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
    numThreads=args.jobs if not (args.geant4 and args.jobs == -1) else 1,
    logLevel=acts.logging.Level(args.loglevel),
    outputDir="" if outputDirFinal is None else str(outputDirFinal),
)

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
            multiplicity=args.gen_multiplicity,
            rnd=rnd,
            # disable output if -o not specified, to match full_chain_itk.py behaviour
            outputDirRoot=outputDirRoot if args.output_dir is not None else None,
            outputDirCsv=outputDirCsv,
        )
    else:
        from acts.examples.simulation import addPythia8

        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=args.gen_nparticles,
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
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
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
        )
    else:
        if s.config.numThreads != 1:
            logger.error(
                f"****** Geant 4 simulation does not support multi-threading (threads={s.config.numThreads}) ******"
            )

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
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
        )

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=args.digi_config,
    rnd=rnd,
    outputDirRoot=outputDirRoot,
    outputDirCsv=outputDirCsv,
)

if args.reco:

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

    if args.itk and args.algorithm == SeedingAlgorithm.Default:
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
        seedingAlgorithm=args.algorithm,
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
        outputDirRoot=outputDirFinalRoot,
        outputDirCsv=outputDirFinalCsv,
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
            outputDirRoot=outputDirFinalRoot,
            outputDirCsv=outputDirFinalCsv,
        )

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
            ) if not args.simple_ckf and args.algorithm != SeedingAlgorithm.TruthSmeared else {}),
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
            ) if args.algorithm != SeedingAlgorithm.TruthSmeared else {}),
            # ITk volumes from Noemi's plot
            pixelVolumes=[8, 9, 10, 13, 14, 15, 16, 18, 19, 20],
            stripVolumes=[22, 23, 24],
            maxPixelHoles=1,
            maxStripHoles=2,
        ) if not args.simple_ckf else CkfConfig()
        # fmt: on

    addCKFTracks(
        s,
        trackingGeometry,
        field,
        trackSelectorConfig=trackSelectorConfig,
        ckfConfig=ckfConfig,
        **(dict(twoWay=False) if not args.simple_ckf else {}),
        **(dict(writeTrackSummary=False) if args.reduce_output >= 1 else {}),
        outputDirRoot=outputDirFinalRoot,
        outputDirCsv=outputDirFinalCsv,
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
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
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
            writeCovMat=True,
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
        )

    elif args.ambi_solver == "greedy":

        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=3,
                maximumIterations=10000,
                nMeasurementsMin=6,
            ),
            **(dict(writeTrackSummary=False) if args.reduce_output >= 1 else {}),
            outputDirRoot=outputDirRoot,
            outputDirCsv=outputDirCsv,
        )

    if args.vertexing:
        addVertexFitting(
            s,
            field,
            vertexFinder=VertexFinder.AMVF,
            outputDirRoot=outputDir,
        )

s.run()
