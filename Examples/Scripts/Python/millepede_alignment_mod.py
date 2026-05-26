#!/usr/bin/env python3

from math import floor
import os
import argparse
import pathlib
import subprocess
from pathlib import Path

import acts
import acts.examples

import enum
import matplotlib.pyplot as plt


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
    discardUnconstrainedTrackPar: bool = True,
    outFileInternalSolving: str = "ActsInternalAlignment_Result.txt",
    milleOutput: str = "MilleBinary.root",
    outFileDecomposition: str = "ActsEigenVecs.txt"
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
        outFileDecomposition=outFileDecomposition
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
    txtOutput="ActsAlignmentRes.txt",
):

    solver = ActsSolverFromMille(
        level=logLevel,
        milleInput=milleInput,
        outFile=txtOutput,
        trackingGeometry=trackingGeometry,
        magneticField=magField,
        fixModules=fixModules,
    )
    s.addAlgorithm(solver)
    return s


# Some helpers for visualising the result
class DoFs(enum.IntEnum):
    dx = 1
    dy = 2
    dz = 3
    rx = 4
    ry = 5
    rz = 6


class AlignmentParResult:
    def __init__(self, theLine: str = " 0 0 0 0 "):
        if len(theLine.split()) == 4:
            labelTxt, ValTxt, seed, n = theLine.split()
            sigmaTxt = " 0. "
        else:
            labelTxt, ValTxt, seed, val2, sigmaTxt, n = theLine.split()
        self.label = int(labelTxt)
        self.value = float(ValTxt)
        self.sigma = float(sigmaTxt)
        self.seed = float(seed)
        self.Counts = int(n)

    def signif(self):
        return (self.value - self.seed) / self.sigma


def readMpEve(fname: str = "millepede.eve"):
    res = {}
    evValue = 0
    with open(fname, "r") as f:
        for line in f.readlines():
            tokens = line.split()
            if len(tokens) == 0:
                continue
            if "Eigenvector" in tokens:
                evIndex = tokens[1]
                evValue = tokens[4]
                res[evValue] = {"EV": evIndex, "Elements": []}
                # print ("EV ",evIndex,"=",evValue)
                continue
            atLabel = True
            atVal = False
            label = 0
            val = 0
            for thing in tokens:
                if atLabel:
                    label = thing
                    atLabel = False
                    atVal = True
                elif atVal:
                    val = thing
                    atVal = False
                    atLabel = True
                else:
                    print("should never get here")
                if label != 0 and atLabel:
                    thisPar = AlignmentParResult()
                    thisPar.label = int(label)
                    thisPar.value = float(val)
                    res[evValue]["Elements"] += [thisPar]
    return res


def readAli(fin):
    content = []
    try: 
        with open(fin) as f:
            for line in f.readlines()[1:]:
                content += [AlignmentParResult(line)]
    except FileNotFoundError: 
        print (f"Could not find '{fin}' - skipping visualisation for a set of results")
    return content


def selByDoF(results: list, DoF: int):
    return [res for res in results if res.label % 6 == int(DoF % 6)]


def plotPar(
    resultFileActs, resultFileDirect, resultFileMP2, DoF, fname, trueValues=[], **kwargs
):

    results_mp2 = readAli(resultFileMP2)
    results_direct = readAli(resultFileDirect)
    results_acts = readAli(resultFileActs)

    res_mp2 = selByDoF(results_mp2, DoF)
    res_direct = selByDoF(results_direct, DoF)
    res_acts = selByDoF(results_acts, DoF)

    fig, ax = plt.subplots()
    labels = [r.label for r in res_mp2]
    values = [r.value for r in res_mp2]
    errorsY = [r.sigma for r in res_mp2]

    direct_labels = [r.label for r in res_direct]
    direct_values = [-1 * r.value for r in res_direct]
    direct_errorsY = [r.sigma for r in res_direct]

    acts_labels = [r.label for r in res_acts]
    # sign convention
    acts_values = [-1 * r.value for r in res_acts]
    acts_errorsY = [r.sigma for r in res_acts]

    ax.errorbar(
        direct_labels,
        direct_values,
        direct_errorsY,
        color="xkcd:light blue",
        marker="s",
        linestyle="None",
        label="ACTS direct",
    )

    ax.errorbar(
        acts_labels,
        acts_values,
        acts_errorsY,
        color="xkcd:orange",
        marker="D",
        linestyle="None",
        fillstyle="none",
        label="ACTS via Mille",
    )
    ax.errorbar(labels, values, errorsY, ActsInternalAlignment_Resultlabel="MillePede-II", **kwargs)

    corrVals = []
    corrLabels = []
    corrErrors = []
    for label, value in trueValues:
        injected = [-1 * value for r in res_mp2 if r.label == label]
        if len(injected) > 0:
            corrVals += injected
            corrLabels += [label]
            corrErrors += [0]
    ax.errorbar(
        corrLabels,
        corrVals,
        corrErrors,
        color="b",
        marker="o",
        fillstyle="none",
        label="injected",
    )
    ax.set_xlabel("Parameter label")
    ax.set_ylabel("Correction")
    ax.grid(which="major", axis="y")
    ax.legend()
    fig.savefig(fname, bbox_inches="tight", pad_inches=0.2)


def plotWM(evFile, fname, evToPlot: list, **kwargs):

    results_EV = readMpEve(evFile)
    cols = ["k", "b", "r", "orange", "green", "violet"]
    fig, ax = plt.subplots(
        2, 3, sharey="row", sharex=True, figsize=[3 * 6.4, 2 * 4.8], dpi=300
    )
    dofs = [DoFs.dx, DoFs.dy, DoFs.dz, DoFs.rx, DoFs.ry, DoFs.rz]
    dofnames = ["tx", "ty", "tz", "rx", "ry", "rz"]
    for i, D in enumerate(dofs):
        ix = i % 3
        iy = floor(i / 3)
        ientry = 0
        for eigenval, dic in results_EV.items():
            ientry += 1
            if ientry - 1 not in evToPlot:
                continue
            theIndex = dic["EV"]
            pars = selByDoF(dic["Elements"], D)

            labels = [r.label for r in pars]
            values = [r.value for r in pars]

            ax[iy, ix].errorbar(
                labels,
                values,
                [0 for v in values],
                label=f"EV {theIndex} = {eigenval}",
                **kwargs,
            )

        ax[iy, ix].set_xlabel("Parameter label")
        ax[iy, ix].set_ylabel(f"Eigenvec {dofnames[i]} component")
        ax[iy, ix].grid(which="major", axis="y")
    ax[1, 2].legend()

    fig.savefig(fname, bbox_inches="tight", pad_inches=0.2)


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
detector = getTelescopeDetector()
trackingGeometry = detector.trackingGeometry()
decorators = detector.contextDecorators()

gcx = acts.GeometryContext.dangerouslyDefaultConstruct()

# inject a known misalignment.

# Misalign the second tracking layer (ID = 4).
layerToBump = acts.GeometryIdentifier(layer=4, volume=1, sensitive=1)

# shift this layer by 200 microns in the global Z direction
leShift = AlignmentGeneratorGlobalShift(shift = acts.Vector3(0, 0, 200.0e-3))

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
# fixModules = {
#     acts.GeometryIdentifier(layer=2, volume=1, sensitive=1),
#     acts.GeometryIdentifier(layer=10, volume=1, sensitive=1),
#     acts.GeometryIdentifier(layer=18, volume=1, sensitive=1),
# }

fixModules = set()

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
# s.addContextDecorator(alignDecoGlobShift)
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
        impactMax=50 * u.mm,
    ),
    # why do we need to specify this here again? Not taken from event context?
    seedFinderOptionsArg=SeedFinderOptionsArg(bFieldInZ=2 * u.T),
    seedingAlgorithm=SeedingAlgorithm.GridTriplet,
    initialSigmas=[
        20 * u.mm,
        20 * u.mm,
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
milleBinary = outputDir / "MyBinary.dat"
# And add our alignment sandbox
actsInternalResFile = "ActsInternalSolver.txt"
addAlignmentSandbox(s, trackingGeometry, field, fixModules, milleOutput=milleBinary,outFileInternalSolving=actsInternalResFile)
actsResultFile = "ActsAlignmentRes.txt"
# And finally read back and solve in ACTS
addSolverFromMille(
    s,
    trackingGeometry,
    field,
    fixModules,
    milleInput=milleBinary,
    txtOutput=actsResultFile,
)

s.run()

### now, we solve the same alignment problem using the Millepede solver

print(" ===== Now aligning with MP2 =====")

# Millepede steering file - this should be sensibly generated in real life use
pedeSteeringArgs = [
    # tell pede to load our binary
    "Cfiles",
    str(milleBinary),
    # set solution method to matrix inversion, max 5 internal iterations, tolerance 0.8
    "Method  diagonalization 5 0.8",
    # "Method  Inversion 5 0.8",
    # print progress
    "monitorprogress 1 10000",
    # min 10 entries per degree of freedom to activate in fit
    "entries 10 10 2",
    # iterations for outlier downweighting
    "outlierdownweighting 4",
    # parallelisation of inversion
    "threads 24",
    # only recalculate global matrix once
    "matiter 1",
    # print statistics
    "printcounts 2",
    # fill built-in residual plots
    "monitorresiduals",
    # count records
    "countrecords",
    # reject tracks with poor condition number
    "maxlocalcond 10.",
    # outlier warning
    "outlierfracwarnthreshold 1",
    # "Parameter",
]
# fix the last layer too
# pedeSteeringArgs+= [f"{i}  0.0   -1.   " for i in range(43,49)]
# lock y-movements
# pedeSteeringArgs+= [f"{i*6+2}  0.0   -1.   " for i in range(8)]
steeringFile = outputDir / "pedeSteerMaster.txt"

with open(steeringFile, "w") as fSteer:
    fSteer.writelines([a + "\n" for a in pedeSteeringArgs])

# now run the fitter
subprocess.run(["pede", str(steeringFile)])

# check the return code
with open("millepede.end", "r") as mpend:
    tokens = mpend.readline().split()
    stat = int(tokens[0])
    descr = " ".join(tokens[1:])

statCodes = {
    (-1, -1): "crash",
    (0, 0): "normal exit",
    (1, 1): "tolerable warnings",
    (2, 5): "serious warnings",
    (10, 40): "abort",
}
res = "unknown"
for interval, name in statCodes.items():
    if interval[0] <= stat and stat <= interval[1]:
        res = name
print(f"Millepede-II Alignment fit ended with {res}.\n   Status code {stat}: {descr}")

# visualise!

plot_raw_dx = plotPar(
    actsResultFile,
    actsInternalResFile,
    "millepede.res",
    DoFs.dx,
    "Results_dx.png",
    # trueValues=[(6 * k+1, 10.e-3) for k in range(9) ],
    fmt=".k",
)
plot_raw_dy = plotPar(
    actsResultFile,
    actsInternalResFile,
    "millepede.res",
    DoFs.dy,
    "Results_dy.png",
    # trueValues=[(6 * k+2, -10.e-3) for k in range(9) ],
    fmt=".k",
)
plot_raw_dz = plotPar(
    actsResultFile,
    actsInternalResFile,
    "millepede.res",
    DoFs.dz,
    "Results_dz.png",
    # trueValues=[(6 * k+3, 20.e-3) for k in range(9) ],
    trueValues=[(3, 0.200)],
    fmt=".k",
)
plot_raw_rx = plotPar(
    actsResultFile,
    actsInternalResFile,
    "millepede.res",
    DoFs.rx,
    "Results_rx.png",
    fmt=".k",
)
plot_raw_ry = plotPar(
    actsResultFile,
    actsInternalResFile,
    "millepede.res",
    DoFs.ry,
    "Results_ry.png",
    fmt=".k",
)
plot_raw_rz = plotPar(
    actsResultFile,
    actsInternalResFile,
    "millepede.res",
    DoFs.rz,
    "Results_rz.png",
    fmt=".k",
)

# evSource = "millepede.eve"
evSource = "ActsEigenVecs.txt"
for EV in range(9):
    plot_ev_single = plotWM(evSource, f"EigenVec_{EV}.png", evToPlot=[EV])

plot_ev_group = plotWM(
    evSource,
    f"EigenVecs_first6.png",
    evToPlot=range(6),
)
for EV in range(9):
    plot_ev_single = plotWM(evSource, f"EigenVec_{EV}.png", evToPlot=[EV])

plot_ev_group = plotWM(
    evSource,
    f"EigenVecs_first4.png",
    evToPlot=range(4),
)
plot_ev_group = plotWM(
    evSource,
    f"EigenVecs_first5.png",
    evToPlot=range(5),
)
plot_ev_group = plotWM(
    evSource,
    f"EigenVecs_first8.png",
    evToPlot=range(8),
)