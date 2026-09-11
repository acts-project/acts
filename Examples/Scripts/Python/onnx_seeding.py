#!/usr/bin/env python3
"""
Full reconstruction chain on the OpenDataDetector using an ONNX-based
PythonCallable seeder. Demonstrates how to plug an ONNX seeding model
(here: GUNTAM) into ACTS via SeedingAlgorithm.PythonCallable.

Usage:
    python onnx_seeding.py --model-path /path/to/model.onnx
    python onnx_seeding.py --model-path /path/to/model.onnx --ttbar --ttbar-pu 200
"""

import argparse
import pathlib

import acts
import acts.examples
from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleGun,
    addPythia8,
    addGenParticleSelection,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    SeedingAlgorithm,
    addSeeding,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory
from guntam_transformer_seeder import guntam_transformer_seeder

u = acts.UnitConstants


parser = argparse.ArgumentParser(
    description="ODD full chain with an ONNX PythonCallable seeder (GUNTAM)"
)
parser.add_argument(
    "--output",
    "-o",
    help="Output directory",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "onnx_seeding_output",
)
parser.add_argument("--events", "-n", help="Number of events", type=int, default=100)
parser.add_argument(
    "--threads",
    help="Number of sequencer threads, i.e. events processed in parallel (-1 = all cores). "
    "For CPU inference ensure --threads * --onnx-threads does not exceed available cores.",
    type=int,
    default=-1,
)
parser.add_argument(
    "--onnx-threads",
    help="ONNX intra- and inter-op thread count for CPU inference (default 1).",
    type=int,
    default=1,
)
parser.add_argument(
    "--gpu",
    help="Use CUDA GPU for ONNX inference (CUDAExecutionProvider with CPU fallback).",
    action="store_true",
)
parser.add_argument(
    "--model-path",
    help="Path to the ONNX seeding model",
    type=pathlib.Path,
    required=True,
)
parser.add_argument(
    "--ttbar",
    help="Use Pythia8 ttbar instead of particle gun",
    action="store_true",
)
parser.add_argument(
    "--ttbar-pu",
    help="Number of pile-up events (only used with --ttbar)",
    type=int,
    default=200,
)
parser.add_argument(
    "--output-root",
    help="Write ROOT output files",
    default=True,
    action=argparse.BooleanOptionalAction,
)
parser.add_argument(
    "--output-csv",
    help="Write CSV output files",
    default=True,
    action=argparse.BooleanOptionalAction,
)

args = parser.parse_args()

outputDir = args.output
outputDir.mkdir(parents=True, exist_ok=True)

geoDir = getOpenDataDetectorDirectory()
actsDir = pathlib.Path(__file__).resolve().parents[3]

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = actsDir / "Examples/Configs/odd-digi-smearing-config.json"
oddSeedingSel = actsDir / "Examples/Configs/odd-seeding-config.json"

oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)
detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
trackingGeometry = detector.trackingGeometry()

field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=args.events,
    numThreads=args.threads,
    outputDir=str(outputDir),
    # ONNX SIMD kernels raise harmless FP underflow signals; disabling avoids false-positive failures unrelated to this script.
    failOnUnmaskedFpe=False,
)

if not args.ttbar:
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-3.0, 3.0),
        PhiConfig(0.0, 360.0 * u.degree),
        ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 1.0 * u.ns),
        ),
        multiplicity=200,
        rnd=rnd,
    )
else:
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=args.ttbar_pu,
        vtxGen=acts.examples.GaussianVertexGenerator(
            mean=acts.Vector4(0, 0, 0, 0),
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
        ),
        rnd=rnd,
        outputDirRoot=outputDir if args.output_root else None,
        outputDirCsv=outputDir if args.output_csv else None,
    )

    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            pt=(150 * u.MeV, None),
        ),
    )

addFatras(
    s,
    trackingGeometry,
    field,
    enableInteractions=True,
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
    rnd=rnd,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=oddDigiConfig,
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        eta=(-3.0, 3.0),
        measurements=(9, None),
        removeNeutral=True,
    ),
)

addSeeding(
    s,
    trackingGeometry,
    field,
    seedingAlgorithm=SeedingAlgorithm.PythonCallable,
    customSeeder=guntam_transformer_seeder,
    customSeederConfig={
        "model_path": str(args.model_path),
        "num_threads": args.onnx_threads,
        "providers": (
            ["CUDAExecutionProvider", "CPUExecutionProvider"] if args.gpu else None
        ),
    },
    geoSelectionConfigFile=oddSeedingSel,
    initialSigmas=[
        1 * u.mm,
        1 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0 * u.e / u.GeV,
        1 * u.ns,
    ],
    initialSigmaQoverPt=0.1 * u.e / u.GeV,
    initialSigmaPtRel=0.1,
    initialVarInflation=[1.0] * 6,
    particleHypothesis=acts.ParticleHypothesis.muon,
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(1.0 * u.GeV if args.ttbar else 0.0, None),
        absEta=(None, 3.0),
        loc0=(-4.0 * u.mm, 4.0 * u.mm),
        nMeasurementsMin=7,
        maxHoles=2,
        maxOutliers=2,
    ),
    CkfConfig(
        chi2CutOffMeasurement=15.0,
        chi2CutOffOutlier=25.0,
        numMeasurementsCutOff=2,
        seedDeduplication=True,
        stayOnSeed=False,
        pixelVolumes=[16, 17, 18],
        stripVolumes=[23, 24, 25],
        maxPixelHoles=1,
        maxStripHoles=2,
        constrainToVolumes=[
            2,  # beam pipe
            32,
            4,  # beam pipe gap
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
    ),
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
    writeCovMat=True,
)

addAmbiguityResolution(
    s,
    AmbiguityResolutionConfig(
        maximumSharedHits=3, maximumIterations=1000000, nMeasurementsMin=7
    ),
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
    writeCovMat=True,
)

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.AMVF,
    outputDirRoot=outputDir if args.output_root else None,
    outputDirCsv=outputDir if args.output_csv else None,
)

s.run()
