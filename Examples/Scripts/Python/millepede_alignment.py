#!/usr/bin/env python3

import os
import argparse
import pathlib
from pathlib import Path

import acts
import acts.examples

from acts.examples import (
    TelescopeDetector,
    AlignedTelescopeDetector,
    Sequencer,
    StructureSelector,
    RandomNumbers,
    GaussianVertexGenerator,
)

from acts.examples.alignment import (
    AlignmentDecorator,
    GeoIdAlignmentStore,
    AlignmentGeneratorGlobalShift,
)
from acts.examples.alignmentmillepede import (
    MillePedeAlignmentSandbox,
    ActsSolverFromMille,
)

from acts.examples.simulation import (
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    ParticleSelectorConfig,
    addParticleGun,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
    SeedingAlgorithm,
    TrackSelectorConfig,
    addSeeding,
    SeedingAlgorithm,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    CkfConfig,
    addCKFTracks,
    TrackSelectorConfig,
)


# Helper to instantiate a telescope detector.
# The square sensors are oriented in the global
# y-direction and cover the x-z plane.
# You can change the number of layers by resizing
# the "bounds", "stereos" and "positions" arrays accordingly.
#
# By default, will be digitised as 25 x 100 pixel grid.
#
# The "layer" field of the geo ID of this detector
# will move in steps of 2 for each module (2,4,6,..,18 for 9 layers).
# Everything is located in volume "1".
#
# In the alignment, at least 4 layers are expected (layer ID 8 will
# be fixed as the alignment reference).
def getTelescopeDetector(misaligned: bool = False):
    bounds = [200, 200]
    positions = [30, 60, 90, 120, 150, 180, 210, 240, 270]
    stereos = [0] * len(positions)
    if not misaligned:
        detector = TelescopeDetector(
            bounds=bounds, positions=positions, stereos=stereos, binValue=1
        )
    else:
        detector = AlignedTelescopeDetector(
            bounds=bounds, positions=positions, stereos=stereos, binValue=1
        )

    return detector


# Add the alignment sandbox algorithm
def addAlignmentSandbox(
    s: Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    magField: acts.MagneticFieldProvider,
    fixModules: set,
    inputMeasurements: str = "measurements",
    inputTracks: str = "ckf_tracks",
    logLevel: acts.logging.Level = acts.logging.INFO,
    milleOutput: str = "MilleBinary.root",
):
    sandbox = MillePedeAlignmentSandbox(
        level=logLevel,
        milleOutput=milleOutput,
        inputMeasurements=inputMeasurements,
        inputTracks=inputTracks,
        trackingGeometry=trackingGeometry,
        magneticField=magField,
        fixModules=fixModules,
    )
    s.addAlgorithm(sandbox)
    return s


def addSolverFromMille(
    s: Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    magField: acts.MagneticFieldProvider,
    fixModules: set,
    logLevel: acts.logging.Level = acts.logging.INFO,
    milleInput: str = "MilleBinary.root",
):

    solver = ActsSolverFromMille(
        level=logLevel,
        milleInput=milleInput,
        trackingGeometry=trackingGeometry,
        magneticField=magField,
        fixModules=fixModules,
    )
    s.addAlgorithm(solver)
    return s


u = acts.UnitConstants

# Can also use zero-field - but not healthy for track covariance.
# field = acts.NullBField()
field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))


parser = argparse.ArgumentParser(
    description="MillePede alignment demo with the Telescope Detector"
)
parser.add_argument(
    "--output",
    "-o",
    help="Output directory",
    type=pathlib.Path,
    default=pathlib.Path.cwd() / "mpali_output",
)
parser.add_argument("--events", "-n", help="Number of events", type=int, default=2000)
parser.add_argument("--skip", "-s", help="Number of events", type=int, default=0)


args = parser.parse_args()

outputDir = args.output
# ensure out output dir exists
os.makedirs(outputDir, exist_ok=True)

# Instantiate the telescope detector - with alignment enabled
detector = getTelescopeDetector(misaligned=True)
trackingGeometry = detector.trackingGeometry()
decorators = detector.contextDecorators()

# inject a known misalignment.

# Misalign the second tracking layer (ID = 4).
layerToBump = acts.GeometryIdentifier(layer=4, volume=1, sensitive=1)

# shift this layer by 200 microns in the global Z direction
leShift = AlignmentGeneratorGlobalShift()
leShift.shift = acts.Vector3(0, 0, 200.0e-3)

# now add some boilerplate code to make this happen
alignDecoConfig = AlignmentDecorator.Config()
alignDecoConfig.nominalStore = GeoIdAlignmentStore(
    StructureSelector(trackingGeometry).selectedTransforms(
        acts.GeometryContext.dangerouslyDefaultConstruct(), layerToBump
    )
)
alignDecoConfig.iovGenerators = [((0, 10000000), leShift)]
alignDeco = AlignmentDecorator(alignDecoConfig, acts.logging.WARNING)
contextDecorators = [alignDeco]

# decide on at least on detector module to fix in place
# as a reference for the alignment.
# By default, fix the innermost layer.
fixModules = {acts.GeometryIdentifier(layer=2, volume=1, sensitive=1)}


# More Boilerplate code - for setting up the sequence

rnd = RandomNumbers(seed=42)
s = Sequencer(
    events=args.events,
    skip=args.skip,
    numThreads=1,
    outputDir=str(outputDir),
)

# Add a context with the alignment shift - sim, digi and
# initial reco will "see" the distorted detector
s.addContextDecorator(alignDeco)

# Run particle gun and fire some muons at our telescope
addParticleGun(
    s,
    MomentumConfig(
        10 * u.GeV,
        100 * u.GeV,
        transverse=True,
    ),
    EtaConfig(-0.3, 0.3),
    # aim roughly along +Y...
    PhiConfig(60 * u.degree, 120 * u.degree),
    ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
    vtxGen=GaussianVertexGenerator(
        mean=acts.Vector4(0, 0, 0, 0),
        stddev=acts.Vector4(5.0 * u.mm, 0.0 * u.mm, 5.0 * u.mm, 0.0 * u.ns),
    ),
    multiplicity=1,
    rnd=rnd,
)
# fast sim
addFatras(
    s,
    trackingGeometry,
    field,
    enableInteractions=True,
    outputDirRoot=None,
    outputDirCsv=None,
    outputDirObj=None,
    rnd=rnd,
)

# digitise with the default 25x100 pixel config
srcdir = Path(__file__).resolve().parent.parent.parent.parent
addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=srcdir / "Examples/Configs/telescope-digi-smearing-config.json",
    outputDirRoot=None,
    outputDirCsv=None,
    rnd=rnd,
)
addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        measurements=(3, None),
        removeNeutral=True,
    ),
)

# Run grid seeder
addSeeding(
    s,
    trackingGeometry,
    field,
    # Settings copied from existing snippet, slightly adapted
    seedFinderConfigArg=SeedFinderConfigArg(
        r=(20 * u.mm, 200 * u.mm),
        deltaR=(1 * u.mm, 300 * u.mm),
        collisionRegion=(-250 * u.mm, 250 * u.mm),
        z=(-100 * u.mm, 100 * u.mm),
        maxSeedsPerSpM=1,
        sigmaScattering=5,
        radLengthPerSeed=0.1,
        minPt=0.5 * u.GeV,
        impactMax=3 * u.mm,
    ),
    # why do we need to specify this here again? Not taken from event context?
    seedFinderOptionsArg=SeedFinderOptionsArg(bFieldInZ=2 * u.T),
    seedingAlgorithm=SeedingAlgorithm.GridTriplet,
    initialSigmas=[
        3 * u.mm,
        3 * u.mm,
        1 * u.degree,
        1 * u.degree,
        0 * u.e / u.GeV,
        1 * u.ns,
    ],
    initialSigmaQoverPt=0.1 * u.e / u.GeV,
    initialSigmaPtRel=0.1,
    initialVarInflation=[1.0] * 6,
    # This file should be adapted if you add layers and want to include
    # them in the seeding.
    geoSelectionConfigFile=srcdir / "Examples/Configs/telescope-seeding-config.json",
    outputDirRoot=None,
)

# Add CKF track finding
addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(),
    CkfConfig(
        chi2CutOffMeasurement=150.0,
        chi2CutOffOutlier=250.0,
        numMeasurementsCutOff=50,
        seedDeduplication=(True),
        stayOnSeed=True,
    ),
    outputDirRoot=None,
    writePerformance=False,
    writeTrackSummary=False,
)

# And add our alignment sandbox
addAlignmentSandbox(
    s, trackingGeometry, field, fixModules, milleOutput=outputDir / "MyBinary.root"
)
# And finally read back and solve in ACTS
addSolverFromMille(
    s, trackingGeometry, field, fixModules, milleInput=outputDir / "MyBinary.root"
)

s.run()
