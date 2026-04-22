# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import argparse

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
        help=("Positional error with path length [%]"),
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
