#!/usr/bin/env python3

# This file is part of the ACTS project.
#
# Copyright (C) 2016 CERN for the benefit of the ACTS project
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.

from pathlib import Path
import argparse
import sys

import acts
import acts.examples
from acts.examples.traccc import *

u = acts.UnitConstants


def make_sequencer(
    s: acts.examples.Sequencer,
    detectorFile: Path,
    digitizationFile: Path,
    bfieldFile: Path,
    dataDirectory: Path,
    conditionsFile: Path,
    materialFile: Path = Path(),
    gridFile: Path = Path(),
):
    """
    Configure the sequencer with the traccc GPU sequence algorithm.

    Parameters
    ----------
    s:
        The sequencer to configure.
    detectorFile:
        Path to the detray detector JSON file.
    digitizationFile:
        Path to the digitization config JSON file.
    bfieldFile:
        Path to the covfie magnetic field binary file.
    dataDirectory:
        Directory containing per-event cell CSV files.
    conditionsFile:
        Path to the detector conditions JSON file.
    materialFile:
        Path to the material map file (optional).
    gridFile:
        Path to the surface grid file (optional).
    """

    s.addAlgorithm(
        TracccSeqAlgorithm(
            level=acts.logging.INFO,
            detectorFile=str(detectorFile),
            digitizationFile=str(digitizationFile),
            bfieldFile=str(bfieldFile),
            dataDirectory=str(dataDirectory),
            conditionsFile=str(conditionsFile) if conditionsFile != Path() else "",
            materialFile=str(materialFile) if materialFile != Path() else "",
            gridFile=str(gridFile) if gridFile != Path() else "",
        )
    )


def main():
    parser = argparse.ArgumentParser(
        description="Run the traccc GPU full-chain sequence over csv cell input files."
    )
    parser.add_argument(
        "--events",
        "-n",
        type=int,
        default=10,
        help="Number of events to process",
    )
    parser.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of events to skip at the start",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path.cwd() / "traccc_output",
        help="Directory for sequencer output",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        required=True,
        help="Directory containing per-event cell CSV files (e.g. event000000000-cells.csv)",
    )
    parser.add_argument(
        "--detector-file",
        type=Path,
        required=True,
        help="Path to detray detector JSON file",
    )
    parser.add_argument(
        "--digitization-file",
        type=Path,
        required=True,
        help="Path to digitization config JSON file",
    )
    parser.add_argument(
        "--bfield-file",
        type=Path,
        required=True,
        help="Path to covfie magnetic field file",
    )
    parser.add_argument(
        "--conditions-file",
        type=Path,
        default=None,
        help="Path to detector conditions file (optional)",
    )
    parser.add_argument(
        "--material-file",
        type=Path,
        default=None,
        help="Path to material map file (optional)",
    )
    parser.add_argument(
        "--grid-file",
        type=Path,
        default=None,
        help="Path to surface grid file (optional)",
    )
    parser.add_argument(
        "--log-level",
        choices=["VERBOSE", "DEBUG", "INFO", "WARNING", "ERROR", "FATAL"],
        default="INFO",
        help="Logging level",
    )
    args = parser.parse_args()

    # Validate required inputs exist
    for label, path in [
        ("detector file", args.detector_file),
        ("digitization file", args.digitization_file),
        ("magnetic field file", args.bfield_file),
        ("conditions file", args.conditions_file),
        ("data directory", args.data_dir),
    ]:
        if not path.exists():
            print(f"ERROR: {label} does not exist: {path}", file=sys.stderr)
            sys.exit(1)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    logLevel = getattr(acts.logging, args.log_level)

    s = acts.examples.Sequencer(
        events=args.events,
        skip=args.skip,
        logLevel=logLevel,
        outputDir=str(args.output_dir),
        trackFpes=False,
    )

    make_sequencer(
        s,
        detectorFile=args.detector_file,
        digitizationFile=args.digitization_file,
        bfieldFile=args.bfield_file,
        dataDirectory=args.data_dir,
        conditionsFile=args.conditions_file or Path(),
        materialFile=args.material_file or Path(),
        gridFile=args.grid_file or Path(),
    )

    s.run()


if __name__ == "__main__":
    main()