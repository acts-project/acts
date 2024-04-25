#!/usr/bin/env python3

import argparse
import pathlib

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)

from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    addCKFTracks,
    TrackSelectorConfig,
)

import acts.examples.reconstruction as reconstruction

from collections import namedtuple
from pathlib import Path

from typing import Optional
from enum import Enum

eta = 4
# eta = 2.5

SeedingAlgorithm = reconstruction.SeedingAlgorithm

DetectorName = Enum("DetectorName", "ODD generic")

HashingMetric = Enum("HashingMetric", "dphi dR")

parser = argparse.ArgumentParser()
parser.add_argument("--mu", type=int, default=0)
parser.add_argument("--bucketSize", type=int, default=0)
parser.add_argument(
    "--nevents",
    type=int,
    # default=1000)
    default=100,
)
parser.add_argument("--maxSeedsPerSpM", type=int, default=1)
parser.add_argument("--seedingAlgorithm", type=str, default="Hashing")
parser.add_argument("--saveFilesSmall", type=bool, default=False)
parser.add_argument("--saveFiles", type=bool, default=False)
parser.add_argument("--AnnoySeed", type=int, default=123456789)
parser.add_argument("--zBins", type=int, default=0)
parser.add_argument("--phiBins", type=int, default=0)
parser.add_argument("--metric", type=str, default="dphi")
args = parser.parse_args()

print(args)

mu = args.mu
bucketSize = args.bucketSize
nevents = args.nevents
saveFiles = args.saveFiles
AnnoySeed = args.AnnoySeed
zBins = args.zBins
phiBins = args.phiBins
maxSeedsPerSpM = args.maxSeedsPerSpM

metric = HashingMetric[args.metric]

seedingAlgorithm = SeedingAlgorithm[args.seedingAlgorithm]


if seedingAlgorithm == SeedingAlgorithm.Default:
    bucketSize = 0
    metric = HashingMetric.dphi
    AnnoySeed = 123456789
    zBins = 0
    phiBins = 0

detector = DetectorName.generic


def extractEnumName(enumvar):
    return str(enumvar).split(".")[-1]


u = acts.UnitConstants


def getActsExamplesDirectory():
    return Path(__file__).parent.parent.parent


Config = namedtuple(
    "Config",
    [
        "mu",
        "bucketSize",
        "maxSeedsPerSpM",
        "detector",
        "seedingAlgorithm",
        "metric",
        "AnnoySeed",
        "zBins",
        "phiBins",
    ],
    defaults=[
        None,
        100,
        1000,
        DetectorName.generic,
        SeedingAlgorithm.Hashing,
        HashingMetric.dphi,
        123456789,
        100_000,
        0,
    ],
)
# https://stackoverflow.com/questions/34269772/type-hints-in-namedtuple
Config.__annotations__ = {
    "mu": int,
    "bucketSize": int,
    "maxSeedsPerSpM": int,
    "detector": DetectorName,
    "seedingAlgorithm": SeedingAlgorithm,
    "metric": str,
    "AnnoySeed": int,
    "zBins": int,
    "phiBins": int,
}

config = Config(
    mu=mu,
    bucketSize=bucketSize,
    maxSeedsPerSpM=maxSeedsPerSpM,
    detector=detector,
    seedingAlgorithm=seedingAlgorithm,
    metric=metric,
    AnnoySeed=AnnoySeed,
    zBins=zBins,
    phiBins=phiBins,
)

actsExamplesDir = getActsExamplesDirectory()

if config.detector == DetectorName.ODD:
    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    geoDir = getOpenDataDetectorDirectory()

    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
    oddSeedingSel = geoDir / "config/odd-seeding-config.json"
    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    detector, trackingGeometry, decorators = getOpenDataDetector(
        odd_dir=geoDir, mdecorator=oddMaterialDeco
    )

    digiConfig = oddDigiConfig

    geoSelectionConfigFile = oddSeedingSel

elif config.detector == DetectorName.generic:
    print("Create detector and tracking geometry")

    detector, trackingGeometry, _ = acts.examples.GenericDetector.create()
    digiConfig = (
        actsExamplesDir
        / "Algorithms/Digitization/share/default-smearing-config-generic.json"
    )
    geoSelectionConfigFile = (
        actsExamplesDir
        / "Algorithms/TrackFinding/share/geoSelection-genericDetector.json"
    )
else:
    exit("Detector not supported")

truthSeedRanges = TruthSeedRanges(
    pt=(1.0 * u.GeV, None), eta=(-eta, eta), nHits=(9, None)
)

CKFptMin = 1.0 * u.GeV

doHashing = config.bucketSize > 0
bucketSize = config.bucketSize
npileup = config.mu
maxSeedsPerSpM = config.maxSeedsPerSpM


def get_dir_config(config: Config):
    global main_dir
    outDir = f"detector_{extractEnumName(config.detector)}"
    outDir += "_output"
    doHashing = config.bucketSize > 0
    if doHashing:
        outDir += "_hashing"

    outDir += f"_mu_{config.mu}"
    if doHashing:
        outDir += f"_bucket_{config.bucketSize}"

    outDir += f"_maxSeedsPerSpM_{config.maxSeedsPerSpM}"

    outDir += f"_seedFinderConfig_{'TrackML'}"

    outDir += f"_seedingAlgorithm_{extractEnumName(config.seedingAlgorithm)}Seeding"
    if doHashing:
        if config.metric != "angular":
            outDir += f"_metric_{extractEnumName(config.metric)}"
        outDir += f"_AnnoySeed_{config.AnnoySeed}"
        if config.zBins != 0:
            outDir += f"_zBins_{config.zBins}"
        if config.phiBins != 0:
            outDir += f"_phiBins_{config.phiBins}"
    return outDir


outDir = get_dir_config(config)

outputDir = pathlib.Path.cwd() / outDir

if not outputDir.exists():
    outputDir.mkdir(parents=True)

config_file = open(outputDir / "config_file.txt", "w")
config_file.write(str(config))
config_file.close()

field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
    events=nevents,
    numThreads=1,
    outputDir=str(outputDir),
    trackFpes=False,
)

addPythia8(
    s,
    hardProcess=["Top:qqbar2ttbar=on"],
    npileup=npileup,
    vtxGen=acts.examples.GaussianVertexGenerator(
        stddev=acts.Vector4(0, 0, 50 * u.mm, 0),
        mean=acts.Vector4(0, 0, 0, 0),
    ),
    rnd=rnd,
    # outputDirRoot=outputDir,
    # outputDirCsv=outputDir if saveFiles else None,
)

addFatras(
    s,
    trackingGeometry,
    field,
    preSelectParticles=ParticleSelectorConfig(
        eta=(-eta, eta), pt=(150 * u.MeV, None), removeNeutral=True
    ),
    enableInteractions=True,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir if saveFiles else None,
    rnd=rnd,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=digiConfig,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir if saveFiles else None,
    rnd=rnd,
)

SeedFinderConfigArg = reconstruction.SeedFinderConfigArg
SeedFinderOptionsArg = reconstruction.SeedFinderOptionsArg
HashingTrainingConfigArg = reconstruction.HashingTrainingConfigArg
HashingAlgorithmConfigArg = reconstruction.HashingAlgorithmConfigArg


import numpy as np

cotThetaMax = 1 / (np.tan(2 * np.arctan(np.exp(-eta))))  # =1/tan(2×atan(e^(-eta)))
seedFinderConfigArg = SeedFinderConfigArg(
    r=(None, 200 * u.mm),  # rMin=default, 33mm
    deltaR=(1 * u.mm, 60 * u.mm),
    collisionRegion=(-250 * u.mm, 250 * u.mm),
    z=(-2000 * u.mm, 2000 * u.mm),
    maxSeedsPerSpM=maxSeedsPerSpM,
    sigmaScattering=5,
    radLengthPerSeed=0.1,
    minPt=500 * u.MeV,
    impactMax=3 * u.mm,
    cotThetaMax=cotThetaMax,  # =1/tan(2×atan(e^(-eta)))
)

seedFinderOptionsArg: SeedFinderOptionsArg = SeedFinderOptionsArg(
    bFieldInZ=1.99724 * u.T
)
hashingTrainingConfigArg: HashingTrainingConfigArg = HashingTrainingConfigArg(
    AnnoySeed=config.AnnoySeed, f=2 if metric == HashingMetric.dR else 1
)
hashingAlgorithmConfigArg: HashingAlgorithmConfigArg = HashingAlgorithmConfigArg(
    bucketSize=config.bucketSize,
    zBins=config.zBins,
    phiBins=config.phiBins,
)
initialSigmas: Optional[list] = [
    1 * u.mm,
    1 * u.mm,
    1 * u.degree,
    1 * u.degree,
    0.1 / u.GeV,
    1 * u.ns,
]

initialVarInflation: Optional[list] = [1.0] * 6

addSeeding(
    s,
    trackingGeometry,
    field,
    seedFinderConfigArg,
    seedFinderOptionsArg,
    hashingTrainingConfigArg,
    hashingAlgorithmConfigArg,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-eta, eta), nHits=(9, None)),
    seedingAlgorithm=seedingAlgorithm,
    geoSelectionConfigFile=geoSelectionConfigFile,
    initialSigmas=initialSigmas,
    initialVarInflation=initialVarInflation,
    outputDirRoot=outputDir,
    outputDirCsv=outputDir if saveFiles else None,
)

logLevel = acts.logging.VERBOSE
rootSpacepointsWriter = acts.examples.RootSpacepointWriter(
    level=logLevel,
    inputSpacepoints="spacepoints",
    filePath=str(outputDir / "spacepoints.root"),
)
s.addWriter(rootSpacepointsWriter)

rootSeedsWriter = acts.examples.RootSeedWriter(
    level=logLevel,
    inputSeeds="seeds",
    filePath=str(outputDir / "seeds.root"),
)
s.addWriter(rootSeedsWriter)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(1.0 * u.GeV, None),
        absEta=(None, eta),
        nMeasurementsMin=6,
    ),
    outputDirRoot=outputDir,
    writeTrajectories=False,
    twoWay=False,
    # outputDirCsv=outputDir,
)

s.run()
