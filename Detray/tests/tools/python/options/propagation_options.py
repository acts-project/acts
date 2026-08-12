# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import argparse

from utils import units

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------

""" Parent parser that contains propagation options """


def propagation_options():

    parser = argparse.ArgumentParser(add_help=False)

    # Navigation options
    parser.add_argument(
        "--min_mask_tol",
        "-min_mtol",
        help=("Min. mask tolerance [mm]"),
        default=1e-05,
        type=float,
    )
    parser.add_argument(
        "--max_mask_tol",
        "-max_mtol",
        help=("Max. mask tolerance [mm]"),
        default=3,
        type=float,
    )
    parser.add_argument(
        "--mask_tol_scalor",
        "-mtol_scalor",
        help=("Scale factor for mask tol."),
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--path_tol", "-ptol", help=("Path tolerance [um]"), default=1, type=float
    )
    parser.add_argument(
        "--overstep_tol",
        "-otol",
        help=("Overstep tolerance [um]"),
        default=-300,
        type=float,
    )
    parser.add_argument(
        "--search_window",
        "-sw",
        nargs=2,
        help=("Surface grid search window."),
        default=[0, 0],
        type=int,
    )
    parser.add_argument(
        "--estimate_scattering_noise",
        "-scatt",
        help=("Open mask tol. die to scattering."),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--n_scattering_stddev",
        "-stddev",
        help=("# standard deviations for scattering noise."),
        default=2,
        type=int,
    )
    parser.add_argument(
        "--accumulated_error",
        "-aerr",
        help=("Positional error with path length [%%]"),
        default=0.0001,
        type=float,
    )

    # Parameter transport options
    parser.add_argument(
        "--min_step_size",
        "-min_step",
        help=("Min. Runge-Kutta step size [mm]"),
        default=0.0001,
        type=float,
    )
    parser.add_argument(
        "--max_step_size",
        "-max_step",
        help=("Max. RKN step size [mm]"),
        default=3.40282e38,
        type=float,
    )
    parser.add_argument(
        "--rk_error_tol",
        "-rk_tol",
        help=("Runge-Kutta tolerance size [mm]"),
        default=0.0001,
        type=float,
    )
    parser.add_argument(
        "--path_limit",
        "-plim",
        help=("Max. path length of a track [m]"),
        default=5,
        type=float,
    )
    parser.add_argument(
        "--bethe_energy_loss",
        "-bethe",
        help=("Use Bethe energy loss"),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--covariance_transport",
        "-cov_trnsp",
        help=("Do covaraiance transport"),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--energy_loss_grad",
        "-egrad",
        help=("Use energy loss gradient"),
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--bfield_grad",
        "-bgrad",
        help=("Use B-field gradient"),
        action="store_true",
        default=False,
    )

    return parser


""" Fill a detray propagation config from the parsed commandline options """


def fill_propagation_config(args, config):

    # Configure the navigator
    navigation = config.navigation
    intersection = navigation.intersection

    if args.search_window is not None:
        if len(args.search_window) != 2:
            raise ValueError(
                "Incorrect surface grid search window. Please provide two integer distances."
            )
        navigation.searchWindow = [args.search_window[0], args.search_window[1]]

    if args.min_mask_tol is not None:
        intersection.minMaskTolerance = args.min_mask_tol * units.mm
    if args.max_mask_tol is not None:
        intersection.maxMaskTolerance = args.max_mask_tol * units.mm
    if args.mask_tol_scalor is not None:
        intersection.maskToleranceScalor = args.mask_tol_scalor
    if args.overstep_tol is not None:
        intersection.overstepTolerance = args.overstep_tol * units.um
    if args.path_tol is not None:
        intersection.pathTolerance = args.path_tol * units.um

    if args.estimate_scattering_noise:
        navigation.estimateScatteringNoise = True

        if args.n_scattering_stddev is not None:
            navigation.nScatteringStddev = args.n_scattering_stddev
        if args.accumulated_error is not None:
            navigation.accumulatedError = args.accumulated_error
    else:
        navigation.estimateScatteringNoise = False

        if args.n_scattering_stddev is not None:
            raise ValueError(
                "Option 'n_scattering_stddev' cannot not be configured unless 'estimate_scattering_noise' is activated"
            )
        if args.accumulated_error is not None:
            raise ValueError(
                "Option 'accumulated_error' cannot not be configured unless 'estimate_scattering_noise' is activated"
            )

    # Configure the stepper
    stepping = config.stepping

    if args.min_step_size is not None:
        stepping.minStepsize = args.min_step_size * units.mm
    if args.max_step_size is not None:
        stepping.stepConstraint = args.max_step_size * units.mm
    if args.rk_error_tol is not None:
        stepping.rkErrorTol = args.rk_error_tol * units.mm
    if args.path_limit is not None:
        stepping.pathLimit = args.path_limit * units.m
    stepping.doCovarianceTransport = args.covariance_transport is True
    if args.bethe_energy_loss is not None:
        stepping.useMeanLoss = args.bethe_energy_loss
    if args.energy_loss_grad is not None:
        stepping.useElossGradient = args.energy_loss_grad
    if args.bfield_grad is not None:
        stepping.useFieldGradient = args.bfield_grad

    return config
