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


""" Fill a uniform track generator config from the parsed arguments """


def fill_track_generator_config(args, config):

    if args.phi_steps is not None:
        config.phiSteps = args.phi_steps
    if args.eta_steps is not None:
        config.etaSteps = args.eta_steps
    if args.random_seed is not None:
        config.seed = args.random_seed
    if args.randomize_charge is not None:
        config.randomizeCharge = args.randomize_charge

    if args.eta_range is not None:
        config.etaRange = (args.eta_range[0], args.eta_range[1])

    if args.pT_range is not None and args.p_range is not None:
        raise ValueError(
            "Transverse and total momentum cannot be specified at the same time"
        )
    if args.pT_range is not None:
        config.pT(args.pT_range[0] * units.GeV)
    elif args.p_range is not None:
        config.pTot(args.p_range[0] * units.GeV)

    return config
