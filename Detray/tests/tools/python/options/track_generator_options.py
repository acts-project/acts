# Detray library, part of the ACTS project (R&D line)
#
# (c) 2024 CERN for the benefit of the ACTS project
#
# Mozilla Public License Version 2.0

import argparse

# ------------------------------------------------------------------------------
# Options parsing
# ------------------------------------------------------------------------------

""" Adds options that are common to all track generators """


def common_track_generator_options(parser):

    parser.add_argument(
        "--random_seed",
        "-seed",
        help=("Seed for the random number generator"),
        default=5489,
        type=int,
    )
    parser.add_argument(
        "--pT_range",
        "-pTr",
        nargs=2,
        help=("Transverse momentum [range] of the test particle [GeV]"),
        type=float,
    )
    parser.add_argument(
        "--p_range",
        "-pr",
        nargs=2,
        help=("Total momentum [range] of the test particle [GeV]"),
        type=float,
    )
    parser.add_argument(
        "--eta_range",
        "-eta",
        nargs=2,
        help=("Eta range of generated tracks"),
        default=[-4, 4],
        type=float,
    )
    parser.add_argument(
        "--randomize_charge",
        "-rand_chrg",
        help=("Randomly flip charge sign per track"),
        action="store_true",
        default=False,
    )

    return parser


""" Parent parser that contains random track generator options """


def random_track_generator_options():

    parser = argparse.ArgumentParser(add_help=False)

    common_track_generator_options(parser)

    parser.add_argument(
        "--n_tracks", "-n", help=("Number of test tracks"), default=100, type=int
    )

    return parser


""" Parent parser that contains uniform track generator options """


def uniform_track_generator_options():

    parser = argparse.ArgumentParser(add_help=False)

    common_track_generator_options(parser)

    # Navigation options
    parser.add_argument(
        "--phi_steps", "-n_phi", help=("Number steps in phi"), default=50, type=int
    )
    parser.add_argument(
        "--eta_steps", "-n_eta", help=("Number steps in eta"), default=50, type=int
    )

    return parser
