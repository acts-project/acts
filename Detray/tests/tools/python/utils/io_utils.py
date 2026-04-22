# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

import json
import os
import sys

""" Read the detector name from geometry json file """


def read_detector_name(geometry_file_name, logging):
    if not os.path.isfile(geometry_file_name):
        logging.error(f"Geometry json file not found! ({geometry_file_name})")
        return "unknown_detector"

    with open(geometry_file_name) as geo_file:
        json_geo = json.loads(geo_file.read())
        det_name = json_geo["header"]["common"]["detector"]

        return det_name


""" Uniform access to the momentum range from the CLI arguments """


def get_p_range(parsed_args, logging):

    if parsed_args.p_range is not None:
        return parsed_args.p_range[0], parsed_args.p_range[1]
    elif parsed_args.pT_range is not None:
        return parsed_args.pT_range[0], parsed_args.pT_range[1]
    else:
        logging.error("No momentum configured")
        sys.exit(1)


""" Add CLI arguments from the detector IO options in args """


def add_detector_io_args(arg_list, parsed_args):

    # Always required
    arg_list.extend(
        [
            "--geometry_file",
            parsed_args.geometry_file,
        ]
    )

    # Optional
    if parsed_args.grid_file:
        arg_list.extend(["--grid_file", parsed_args.grid_file])

    if parsed_args.material_file:
        arg_list.extend(["--material_file", parsed_args.material_file])


""" Add CLI arguments from the track generator options in args """


def add_track_generator_args(arg_list, parsed_args):

    arg_list.extend(
        [
            "--random_seed",
            str(parsed_args.random_seed),
            "--eta_range",
            str(parsed_args.eta_range[0]),
            str(parsed_args.eta_range[1]),
        ]
    )

    has_p_range = parsed_args.p_range is not None
    has_pT_range = parsed_args.pT_range is not None
    if has_p_range and has_pT_range:
        print(
            "ERROR: Cannot set total momentum and transverse momentum at the same time"
        )
        sys.exit(1)
    elif has_p_range:
        arg_list.extend(
            ["--p_range", str(parsed_args.p_range[0]), str(parsed_args.p_range[1])]
        )
    elif has_pT_range:
        arg_list.extend(
            ["--pT_range", str(parsed_args.pT_range[0]), str(parsed_args.pT_range[1])]
        )
    else:
        # Random track generator
        arg_list.extend(["--p_range", "1", "1"])

    if parsed_args.randomize_charge:
        arg_list.extend(["--randomize_charge"])

    if hasattr(parsed_args, "n_tracks"):
        # Random track generator
        arg_list.extend(["--n_tracks", str(parsed_args.n_tracks)])
    else:
        # Uniform track generator
        arg_list.extend(
            [
                "--phi_steps",
                str(parsed_args.phi_steps),
                "--eta_steps",
                str(parsed_args.eta_steps),
            ]
        )


""" Add CLI arguments from the propagation options in args """


def add_propagation_args(arg_list, parsed_args):

    arg_list.extend(
        [
            "--min_mask_tolerance",
            str(parsed_args.min_mask_tol),
            "--max_mask_tolerance",
            str(parsed_args.max_mask_tol),
            "--mask_tolerance_scalor",
            str(parsed_args.mask_tol_scalor),
            "--overstep_tolerance",
            str(parsed_args.overstep_tol),
            "--path_tolerance",
            str(parsed_args.path_tol),
            "--search_window",
            str(parsed_args.search_window[0]),
            str(parsed_args.search_window[1]),
            "--rk-tolerance",
            str(parsed_args.rk_error_tol),
            "--path_limit",
            str(parsed_args.path_limit),
            "--minimum_stepsize",
            str(parsed_args.min_step_size),
            "--step_contraint",
            str(parsed_args.max_step_size),
        ]
    )

    if parsed_args.covariance_transport:
        arg_list.extend(["--covariance_transport"])

    if parsed_args.bethe_energy_loss:
        arg_list.extend(["--mean_energy_loss"])

    if parsed_args.energy_loss_grad:
        arg_list.extend(["--eloss_gradient"])

    if parsed_args.bfield_grad:
        arg_list.extend(["--bfield_gradient"])

    if parsed_args.estimate_scattering_noise:
        arg_list.extend(["--estimate_scattering_noise"])
        arg_list.extend(
            [
                "--n_scattering_stddev",
                str(parsed_args.n_scattering_stddev),
                "--accumulated_error",
                str(parsed_args.accumulated_error),
            ]
        )
