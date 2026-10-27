#!/usr/bin/env python3

#########################
# Demo of how to determine and
# apply equality constraints when
# running alignment fits with
# Millepede-II.
#########################

import os
import argparse
from pathlib import Path

import acts
import acts.mille
import acts.examples

from alignment_visualisation import plotAlignImpacts
import numpy as np

from acts.examples import (
    Sequencer,
)
from acts.mille import (
    readMillePedeResult,
    runMillePedeSolver,
    generateMillePedeSteeringFile,
    MillePedeEqualityConstraint,
    MillePedeParameterSeed,
    MillePedeSolutionStrategy,
    MillePedeSteeringConfig,
    MillePedeSolverConfig,
)

from millepede_alignment import (
    addAlignmentSandbox,
    addMisalignmentDeco,
    addTrackingSequence,
    getTelescopeDetector,
)

from telescope_alignment_visualiser import plot_all_modes


# Helper class to express an eigenvector in the Millepede alignment
class mpEigenVec:
    def __init__(
        self,
        index: int,
        EV: float,
        labels=None,
        components=None,
        component_dic: dict = {},
    ):
        self.index = index
        self.EV = EV
        if labels is not None and components is not None:
            self.labels = labels
            self.components = components
        elif len(component_dic) > 0:
            self.labels = np.array([l for l in component_dic.keys()], dtype=int)
            self.components = np.array([v for v in component_dic.values()], dtype=float)

    # for printout purposes, suppress all minor components so that the dominant direction
    # becomes immediately visible to the reader
    def prune(self, threshold=0.01):
        relevant = np.where(np.abs(self.components) > threshold)
        return mpEigenVec(
            self.index,
            self.EV,
            self.labels[relevant],
            self.components[relevant]
            / np.linalg.vector_norm(self.components[relevant]),
        )

    # print helper
    def print(self):
        print(
            f"EV # {self.index} with EV {self.EV:.2f} has labels \n   {self.labels}\n and components \n   {np.round(self.components,3)}\n "
        )


# helper function - read back eigenvectors written by MP-II
def readEvFile(fname):
    # The format consists of a header "Eigenvector <INDEX> with eigenvalue <value> "
    # followed by a set of "label  value" pairs, multiple possible pert line.
    with open(fname, "r") as fin:
        evs = []
        theIndex = None
        theVal = None
        elements = {}
        for li in fin.readlines():
            spl = li.split()
            # header line
            if "Eigenvector" in li:
                if theIndex is not None and theVal is not None:
                    evs.append(mpEigenVec(theIndex, theVal, component_dic=elements))
                    elements.clear()
                theIndex = int(spl[1])
                theVal = float(spl[4])
            # part of an ongoing eigenvector - account for multiple pairs per line
            elif len(spl) > 1:
                elements.update(
                    {
                        int(spl[2 * k]): float(spl[2 * k + 1])
                        for k in range(int(len(spl) / 2))
                    }
                )
    return evs


# package results for visualisation purposes and call
# the dedicated plotter
def drawResult(
    paramsMP2,
    constraints,
    truthInjected: dict,
    fname="AlignmentResults",
    title=None,
    xLim=None,
):
    # package our dictionary of injections to a vector
    injectedVals = np.zeros(54, dtype=float)
    for label, value in truthInjected.items():
        injectedVals[label - 1] = value

    if len(constraints) > 0:
        # translate the injected result
        # into the sub-space permitted
        # by the constraints
        blindSpace = np.zeros(shape=(54, len(constraints)), dtype=float)
        # build the matrix representing the blind
        # directions (constraints as columns)
        for i, constraint in enumerate(constraints):
            for label, value in constraint.labelsAndWeights:
                blindSpace[label - 1, i] = value
        # and then perform the projection
        injectedVals = injectedVals - blindSpace @ np.linalg.solve(
            blindSpace.T @ blindSpace, blindSpace.T @ injectedVals
        )
    # read the fit results into a 3D array
    mp2Info = np.array([[p.label, p.val, p.sigma] for p in paramsMP2], dtype=float)
    # sort ascending by labels
    mp2Info = mp2Info[np.argsort(mp2Info[:, 0])]
    # and extract the columns with the labels / values / uncertainties.
    # For each of those, make one 2D array holding (fit, injected)
    labels = mp2Info[:, 0]
    central = np.zeros(shape=(2, len(labels)), dtype=float)
    sigma = np.zeros(shape=(2, len(labels)), dtype=float)

    # column 0 = fit result
    central[0] = mp2Info[:, 1] * -1  # sign convention different between MP2 and ACTS

    # fill column 1 with the injection
    if truthInjected is not None:
        central[1, :] = injectedVals

    # also update uncertainties, if we have any. Only for fit result
    sigmaToUpdate = np.where(np.isfinite(mp2Info[:, 2])) and np.where(mp2Info[:, 2] > 0)
    sigma[0, sigmaToUpdate] = mp2Info[sigmaToUpdate, 2]

    # lazily hardcode an "okay-looking" format
    formatDict = [
        {
            "label": "Millepede-II",
            "fmt": "o",
            "color": "black",
            "markersize": 5,
        },
        {
            "label": "Injected" + ", Projected" if len(constraints) > 0 else "",
            "marker": "s",
            "markeredgecolor": "#f31d1d",
            "markerfacecolor": "none",
            "linestyle": "none",
            "markersize": 11,
            "markeredgewidth": 1.6,
        },
    ]

    # and hand over to the plotter
    plotAlignImpacts(
        labels,
        central,
        sigma,
        formatDict,
        fname=fname,
        title=title,
        xLim=xLim,
    )


def main():

    logger = acts.getDefaultLogger("MillePedeResultReader", acts.logging.INFO)
    parser = argparse.ArgumentParser(
        description="MillePede alignment demo with the Telescope Detector and equality constraints"
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output directory",
        type=Path,
        default=Path.cwd() / "mpali_output",
    )
    parser.add_argument(
        "--events", "-n", help="Number of events", type=int, default=2000
    )
    parser.add_argument("--skip", "-s", help="Number of events", type=int, default=0)

    args = parser.parse_args()

    # book an output folder
    outputDir = args.output
    os.makedirs(outputDir, exist_ok=True)

    # Set file locations for results
    milleBinaryName = outputDir / "MilleBinary.root"
    millePedeSolution = outputDir / "MillePedeAlignment.txt"
    millePedeSolutionSoftMode = outputDir / "MillePedeAlignment.SoftModeCut.txt"
    millePedeEigenvectors = outputDir / "MillePedeEV.txt"
    # More Boilerplate code - for setting up the sequence

    s = Sequencer(
        events=args.events,
        skip=args.skip,
        numThreads=8,
        outputDir=str(outputDir),
    )

    # Can also use zero-field - but not healthy for track covariance.
    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * acts.UnitConstants.T))

    # build the detector
    detector = getTelescopeDetector()
    trackingGeometry = detector.trackingGeometry()

    # Add a context with the alignment shift - sim and digi will "see" the
    # distorted detector, reconstruction will not
    addMisalignmentDeco(detector, trackingGeometry, s)

    # Add single particle generation, FATRAS, and reco
    # via CKF
    addTrackingSequence(s, field, detector)

    # And add our alignment sandbox
    addAlignmentSandbox(
        s,
        trackingGeometry,
        field,
        fixModules=set(),
        milleOutput=str(milleBinaryName),
    )

    # We do not schedule the solver algorithm.
    # Instead, let's use the MP tooling directly
    # in python via the bindings!

    # run the sequence - generates alignment binary
    s.run()

    ##### First solver run.
    # We pretend we don't know what our tracker is doing.
    # In such a case, it can be instructive to run an
    # Eigendecomposition of the alignment matrix to
    # probe for singular and weak modes.
    # We will do so and then constrain them in a
    # second fit.
    ###########################
    print("\n=== Running initial fit to discover singular modes === ")

    # We already know a trivial 3D translation
    # or rotation of the entire tracker is singular.
    # To save time, lock it down already in
    # iteration 1.
    constraintsStep1 = []
    for j in range(6):
        # build constraints: The sum of all layer movements
        # along each of the three axes must be zero
        # (equivalent: centre of gravity is fixed).
        eqc = MillePedeEqualityConstraint()
        eqc.labelsAndWeights = [(1 + 6 * k + j, 1) for k in range(9)]
        eqc.constraint = 0.0
        constraintsStep1.append(eqc)

    # configure the solver and steering file for Millepede.
    # For this example, we use the defaults for everything and
    # only touch the file names
    steerCfgStep1 = MillePedeSteeringConfig(
        inputFiles=[str(milleBinaryName)],
        # the "Diagonalization" strategy performs eigenvector decomposition,
        # which is the diagnostic tool we want.
        strategy=MillePedeSolutionStrategy.Diagonalization,
        # tell MP to give up after one iteration - we do not yet expect
        # it to be able to converge given singular modes.
        extraLines=["subido"],
        constraints=constraintsStep1,
    )
    # generate the steering file
    fnameSteerStep1 = "steerStep1.txt"
    generateMillePedeSteeringFile(fnameSteerStep1, steerCfgStep1, logger)

    # and configure the solver to use our steering.
    # Also tell it that it should relocate
    # the eigenvector file to a place where
    # we can pick it up for analysis.
    solverCfgStep1 = MillePedeSolverConfig(
        steeringFile=fnameSteerStep1,
        workDir="mpTmpConstraints/Step1",
        evFileName=str(millePedeEigenvectors),
        redirectStdout=str(outputDir / "MillepedeAlignmentStep1.log"),
    )
    # and run the first alignment fit
    runMillePedeSolver(solverCfgStep1, logger)

    # now read back the eigenvectors
    theEVset = readEvFile(millePedeEigenvectors)

    # print the least constrained modes - in a
    # compressed form limited to the leading components.
    # The cut-off value to use can be adjusted - in this example,
    # all "good" eigenmodes are well above 500, while
    # the three remaining singular modes are all below 1.
    theFilteredEv = [e for e in theEVset if e.EV < 500.0]
    print("\n=== Found the following weak / singular modes: ===")
    for thing in theFilteredEv:
        thing.prune(0.03).print()

    # and also make plots of the least constrained eigenmodes
    plot_all_modes(theFilteredEv)
    print("Please see also 'alignment_mode_XYZ.png' in your working dir")

    # Here, we would usually analyse the result and try to understand the Eigenmodes.
    # In this example, we pretend to have already done this.
    # If you look at the printouts or plots, you will see
    # three additional near-singular modes in addition to
    # the three translations we manually vetoed:
    # - a longitudinal telescope mode along y
    # - two transverse telescope (or shear) modes along x/z as function of y

    # For fun, you can disable the constraints in step 1 and see
    # the global translation / rotation modes appear!

    # now form constraints from the additional discovered modes.
    constraints = []
    for EV in theFilteredEv:
        constraints.append(MillePedeEqualityConstraint())
        # this way, we can automate the constraint creation
        # using the prior problem analysis.
        constraints[-1].labelsAndWeights = [
            (EV.labels[i], EV.components[i]) for i in range(len(EV.components))
        ]
        constraints[-1].constraint = 0.0

    # Now we can run the second and final alignment fit.

    print("\n=== Running second fit to perform the alignment with constraints === ")
    # New steering config - add the newly found extra constraints.
    # This time, we can use the standard "inversion" solver.
    steerCfgStep2 = MillePedeSteeringConfig(
        inputFiles=[str(milleBinaryName)],
        strategy=MillePedeSolutionStrategy.Inversion,
        constraints=constraints + constraintsStep1,
    )
    # generate the steering file
    fnameSteerStep2 = "steerStep2.txt"
    generateMillePedeSteeringFile(fnameSteerStep2, steerCfgStep2, logger)

    # and configure the solver
    solverCfgStep2 = MillePedeSolverConfig(
        steeringFile=fnameSteerStep2,
        workDir="mpTmpConstraints/Step2",
        resFileName=str(millePedeSolution),
        redirectStdout=str(outputDir / "MillepedeAlignmentStep2.log"),
    )
    # run the final fit
    runMillePedeSolver(solverCfgStep2, logger)

    ###########
    # Finally, let's try an alternative approach:
    # We can effectively suppress singular / weak
    # modes by adding a penalty term to the objective
    # function that inhibits excessive movement.
    # ATLAS calls this soft mode cut, Millepede
    # has a "regularization"
    # setting that is a close analogue.

    print("\n=== Running alignment with soft mode cut === ")

    # Seed all fit parameters with a small pre-sigma that will
    # constrain them. Use different values for translations
    # and rotations to account for their different units
    preSigmas = []
    for label in range(54):
        s = MillePedeParameterSeed()
        s.label = 1 + label
        # 100 microns or 2 mrad
        s.preSigma = 0.1 if label % 6 < 3 else 0.002
        preSigmas.append(s)

    # Third steering config - this time using pre-sigmas
    # for regularisation.
    # We again use the standard "inversion" solver.
    steerCfgSoftModeCut = MillePedeSteeringConfig(
        inputFiles=[str(milleBinaryName)],
        strategy=MillePedeSolutionStrategy.Inversion,
        # this enables regularisation.
        # First parameter: Strength
        # Second parameter: Default pre-sigma.
        # The penalty term will take the form
        # <strength> x sum (parameter^2)/pre_sigma^2
        # Here, we use the per-parameter pre-sigmas
        # from above rather than the default value.
        extraLines=["regularisation 10  1."],
        seeds=preSigmas,
    )
    # generate the steering file
    fnameSteerSoftModeCut = "steerSoftModeCut.txt"
    generateMillePedeSteeringFile(fnameSteerSoftModeCut, steerCfgSoftModeCut, logger)

    # and configure the solver
    solverCfgSoftMode = MillePedeSolverConfig(
        steeringFile=fnameSteerSoftModeCut,
        workDir="mpTmpConstraints/SoftModeCut",
        resFileName=str(millePedeSolutionSoftMode),
        redirectStdout=str(outputDir / "MillepedeAlignmentSoftMode.log"),
    )
    # run the final fit
    runMillePedeSolver(solverCfgSoftMode, logger)

    # read the results - all solvers follow the Millepede convention,
    # so we can read them the same way
    res_MPII_constraints = readMillePedeResult(str(millePedeSolution), logger)
    res_MPII_softMode = readMillePedeResult(str(millePedeSolutionSoftMode), logger)

    # visualise the outcome!
    drawResult(
        res_MPII_constraints,
        constraints + constraintsStep1,
        {9: 0.05},
        xLim=(-0.2, 0.2),
        fname="AlignmentWithConstraints",
    )
    drawResult(
        res_MPII_softMode,
        constraintsStep1,
        {9: 0.05},
        xLim=(-0.2, 0.2),
        fname="AlignmentWithSoftMode",
    )


if __name__ == "__main__":
    main()
