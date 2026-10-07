#!/usr/bin/env python3

from typing import Any
import os
import argparse
import pathlib
from pathlib import Path

import acts
import acts.mille
import acts.examples

from alignment_visualisation import plotAlignImpacts
import numpy as np

from acts.examples import (
    TelescopeDetector,
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
    MillePedeSolverAlgorithm,
)
from acts.mille import (
    MillePedeSolver,
    readMillePedeResult,
    MillePedeParameterResult,
    generateMillePedeSteeringFile,
    MillePedeEqualityConstraint,
    MillePedeSolutionStrategy,
    MillePedeSteeringConfig,
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
def getTelescopeDetector():
    bounds = [200, 200]
    positions = [30, 60, 90, 120, 150, 180, 210, 240, 270]
    stereos = [0] * len(positions)
    detector = TelescopeDetector(
        bounds=bounds, positions=positions, stereos=stereos, rotDirection=1
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
    discardUnconstrainedTrackPar: bool = True,
    outFileInternalSolving: str = "ActsInternalAlignment_Result.txt",
    outFileDecomposition: str = "ActsInternalAlignment_Eigenvals.txt",
):
    sandbox = MillePedeAlignmentSandbox(
        level=logLevel,
        milleOutput=milleOutput,
        inputMeasurements=inputMeasurements,
        inputTracks=inputTracks,
        trackingGeometry=trackingGeometry,
        magneticField=magField,
        fixModules=fixModules,
        discardUnconstrainedTrackPar=discardUnconstrainedTrackPar,
        outFileInternalSolving=outFileInternalSolving,
        outFileDecomposition=outFileDecomposition,
    )
    s.addAlgorithm(sandbox)
    return s


def addMillePedeSolver(
    s: Sequencer,
    steerCfg: MillePedeSteeringConfig,
    solverCfg: MillePedeSolver.Config,
    logLevel: acts.logging.Level = acts.logging.INFO,
):
    solver = MillePedeSolverAlgorithm(
        level=logLevel, solverConfig=solverCfg, steeringConfig=steerCfg
    )
    s.addAlgorithm(solver)
    return s


def addSolverFromMille(
    s: Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    magField: acts.MagneticFieldProvider,
    fixModules: set,
    logLevel: acts.logging.Level = acts.logging.INFO,
    milleInput: str = "MilleBinary.root",
    outFile: str = "ActsAlignmentViaMille.txt",
):

    solver = ActsSolverFromMille(
        level=logLevel,
        milleInput=milleInput,
        trackingGeometry=trackingGeometry,
        magneticField=magField,
        fixModules=fixModules,
        outFile=outFile,
    )
    s.addAlgorithm(solver)
    return s


def addMisalignmentDeco(detector, trackingGeometry, seq):
    # Instantiate the telescope detector - with alignment enabled
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
    alignDecoConfig.target = AlignmentDecorator.Target.eSim
    alignDeco = AlignmentDecorator(alignDecoConfig, acts.logging.WARNING)

    seq.addContextDecorator(alignDeco)


def drawResult(
    paramsActs,
    paramsActsMille,
    paramsMP2,
    truthInjected=None,
    fname="AlignmentResults",
    title=None,
    xLim=None,
):
    assert len(paramsActs) == len(paramsActsMille) and len(paramsActs) == len(paramsMP2)

    # sort by label index
    actsInfo = np.array([[p.label, p.val, p.sigma] for p in paramsActs], dtype=float)
    actsMilleInfo = np.array(
        [[p.label, p.val, p.sigma] for p in paramsActsMille], dtype=float
    )
    mp2Info = np.array([[p.label, p.val, p.sigma] for p in paramsMP2], dtype=float)

    actsInfo = actsInfo[np.argsort(actsInfo[:, 0])]
    actsMilleInfo = actsMilleInfo[np.argsort(actsMilleInfo[:, 0])]
    mp2Info = mp2Info[np.argsort(mp2Info[:, 0])]
    labels = actsInfo[:, 0]
    central = np.zeros(shape=(4, len(labels)), dtype=float)
    sigma = np.zeros(shape=(4, len(labels)), dtype=float)
    central[0] = actsInfo[:, 1]
    central[1] = actsMilleInfo[:, 1]
    central[3] = mp2Info[:, 1] * -1  # sign convention different
    if truthInjected is not None:
        for k, v in truthInjected.items():
            central[2, np.where(labels == k)] = v

    sigmaToUpdate = np.where(np.isfinite(actsInfo[:, 2])) and np.where(
        actsInfo[:, 2] > 0
    )
    sigma[0, sigmaToUpdate] = actsInfo[sigmaToUpdate, 2]

    sigmaToUpdate = np.where(np.isfinite(actsMilleInfo[:, 2])) and np.where(
        actsMilleInfo[:, 2] > 0
    )
    sigma[1, sigmaToUpdate] = actsMilleInfo[sigmaToUpdate, 2]

    sigmaToUpdate = np.where(np.isfinite(mp2Info[:, 2])) and np.where(mp2Info[:, 2] > 0)
    sigma[3, sigmaToUpdate] = mp2Info[sigmaToUpdate, 2]
    formatDict = [
        {"label": "ACTS-internal", "fmt": "o", "color": "k", "markersize": 5},
        {"label": "ACTS via Mille", "fmt": "o", "color": "blue", "markersize": 5},
        {
            "label": "Injected",
            "marker": "s",
            "markeredgecolor": "#f31d1d",
            "markerfacecolor": "none",
            "linestyle": "none",
            "markersize": 11,
            "markeredgewidth": 1.6,
        },
        {"label": "Millepede-II", "fmt": "o", "color": "red", "markersize": 5},
    ]

    plotAlignImpacts(labels, central, sigma, formatDict, fname, title, xLim=xLim)


def main():
    logger = acts.getDefaultLogger("MillePedeResultReader", acts.logging.INFO)
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
    parser.add_argument(
        "--events", "-n", help="Number of events", type=int, default=2000
    )
    parser.add_argument("--skip", "-s", help="Number of events", type=int, default=0)

    args = parser.parse_args()

    outputDir = args.output
    # ensure out output dir exists
    os.makedirs(outputDir, exist_ok=True)

    # decide on at least on detector module to fix in place
    # as a reference for the alignment.
    # By default, fix the innermost layer.
    fixModules = {
        acts.GeometryIdentifier(layer=2, volume=1, sensitive=1),
        acts.GeometryIdentifier(layer=10, volume=1, sensitive=1),
        acts.GeometryIdentifier(layer=18, volume=1, sensitive=1),
    }

    # Set file locations for results

    milleBinaryName = outputDir / "MilleBinary.root"
    actsInternalSolution = outputDir / "ActsInternalAlignment.txt"
    actsViaMilleSolution = outputDir / "ActsViaMilleAlignment.txt"
    millePedeSolution = outputDir / "MillePedeAlignment.txt"
    # configure the solver and steering file for Millepede.
    # For this example, we use the defaults for everything and
    # only touch the file names
    solverCfg = MillePedeSolver.Config(
        steeringFile="mpsteer.txt",
        workDir="mpTmp",
        resFileName=str(millePedeSolution),
        redirectStdout=str(outputDir / "MillepedeAlignment.log"),
    )

    # EXAMPLE: This is how you could introduce global
    # equality constraints to prevent a 3D
    # translation of the entire detector:
    # The sum of all global movements
    # (labels = 6 * sensor_index
    #           + 0...6 (counts degrees of freedom)
    #           + 1 )
    # is required to be 0 in each
    # translation directions (1..3)
    # Will not apply this as it would be inconsistent
    # with the ACTS solver setup, and we want to
    # compare side-by-side.

    # eqc = [
    #     MillePedeEqualityConstraint(
    #     labelsAndWeights = [(1 + 6 * k+ j,1) for k in range(9) ],
    #     constraint = 0.
    # ) for j in range(3)]

    steerCfg = MillePedeSteeringConfig(
        inputFiles=[str(milleBinaryName)],
        # constraints = eqc             # this would enable constraints
    )

    # More Boilerplate code - for setting up the sequence

    rnd = RandomNumbers(seed=42)
    s = Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=8,
        outputDir=str(outputDir),
    )

    detector = getTelescopeDetector()
    trackingGeometry = detector.trackingGeometry()

    # Add a context with the alignment shift - sim and digi will "see" the
    # distorted detector, reconstruction will not
    addMisalignmentDeco(detector, trackingGeometry, s)

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
        geoSelectionConfigFile=srcdir
        / "Examples/Configs/telescope-seeding-config.json",
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
        s,
        trackingGeometry,
        field,
        fixModules,
        milleOutput=str(milleBinaryName),
        outFileInternalSolving=str(actsInternalSolution),
    )
    # And finally read back and solve in ACTS
    addSolverFromMille(
        s,
        trackingGeometry,
        field,
        fixModules,
        milleInput=str(milleBinaryName),
        outFile=str(actsViaMilleSolution),
    )

    addMillePedeSolver(
        s,
        solverCfg=solverCfg,
        steerCfg=steerCfg,
    )

    s.run()

    # read the results - all solvers follow the Millepede convention,
    # so we can read them the same way
    res_actsInternal = readMillePedeResult(str(actsInternalSolution), logger)
    res_actsViaMille = readMillePedeResult(str(actsViaMilleSolution), logger)
    res_MillepedeII = readMillePedeResult(str(millePedeSolution), logger)

    # visualise the outcome!
    drawResult(res_actsInternal, res_actsViaMille, res_MillepedeII, {3: 0.20})


if __name__ == "__main__":
    main()
